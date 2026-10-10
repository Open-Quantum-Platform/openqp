import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import Mock
import sys
sys.path.insert(0, str(Path(__file__).parent))
from openqp_sync_policy import sync_gateway, INITIAL_SHA, INITIAL_BRANCH, INITIAL_GATEWAY

spec = importlib.util.spec_from_file_location('sync', Path(__file__).with_name('sync.py'))
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)

class Gates(unittest.TestCase):
    def test_review_reuse_rejects_spoofed_or_unrelated_mrs(self):
        mr={'iid':7, 'source_project_id':19,'target_project_id':19,'target_branch':'main',
            'author':{'id':64},'source_branch':'gateway-sync/'+('a'*40)+'/'+('b'*40)}
        self.assertEqual(sync_gateway(19,mr),'a'*40)
        for key,value in [('author',{'id':4}),('source_project_id',17),
                          ('target_project_id',17),('target_branch','feature'),
                          ('source_branch','gateway-sync/fake')]:
            bad=dict(mr);bad[key]=value
            self.assertIsNone(sync_gateway(19,bad))
        self.assertIsNone(sync_gateway(18,mr))

    def test_initial_exception_requires_pinned_mr_author_branch_and_sha(self):
        mr={'iid':4,'source_project_id':19,'target_project_id':19,'target_branch':'main',
            'author':{'id':61},'source_branch':INITIAL_BRANCH,'sha':INITIAL_SHA}
        self.assertEqual(sync_gateway(19,mr),INITIAL_GATEWAY)
        mr['sha']='a'*40
        self.assertIsNone(sync_gateway(19,mr))

    def test_ci_requires_all_platform_jobs_and_exact_private_mr_pipeline(self):
        sha='a'*40
        mr={'sha':sha,'head_pipeline':{'sha':sha,'project_id':19,'source':'merge_request_event','status':'success'}}
        jobs=[{'name':n,'status':'success'} for n in s.REQUIRED_JOBS]
        self.assertEqual(s.pipeline_gate(mr,jobs),'clean')
        self.assertEqual(s.pipeline_gate(mr,jobs[:-1]),'blocked_missing_ci_jobs')
        jobs[0]['status']='skipped'
        self.assertEqual(s.pipeline_gate(mr,jobs),'blocked_incomplete_ci_jobs')
        for key,value in [('sha','b'*40),('project_id',17),('source','push'),('status','failed')]:
            old=mr['head_pipeline'][key];mr['head_pipeline'][key]=value
            self.assertEqual(s.pipeline_gate(mr,jobs),'waiting_successful_ci')
            mr['head_pipeline'][key]=old

    def test_outbound_write_is_rejected_before_credentials_are_read(self):
        with tempfile.TemporaryDirectory() as d:
            obj=s.Sync(Path(d))
            with self.assertRaises(ValueError):obj.api('POST','/projects/17/merge_requests',{})
            with self.assertRaises(ValueError):obj.push_branch('a'*40,'main')
            with self.assertRaises(ValueError):obj.push_branch('a'*40,'unrelated-feature')

class MRGuards(unittest.TestCase):
    def test_merge_commit_policy_is_rejected_before_fetching(self):
        with tempfile.TemporaryDirectory() as d:
            obj=s.Sync(Path(d))
            obj.rotate_token=Mock()
            obj.api=Mock(return_value={'path_with_namespace':s.PATHS[19],
                'merge_method':'merge', 'squash_option':'default_off',
                'only_allow_merge_if_pipeline_succeeds':True,
                'only_allow_merge_if_all_discussions_are_resolved':True})
            obj.fetch=Mock(side_effect=AssertionError('Rejected policy must stop first'))
            self.assertEqual(obj.run()['status'],'blocked_project_policy')

    def test_verified_sync_does_not_request_duplicate_review(self):
        with tempfile.TemporaryDirectory() as d:
            obj=s.Sync(Path(d), dry_run=True)
            sha='a'*40
            mr={'iid':7,'source_project_id':19,'target_project_id':19,
                'target_branch':'main','source_branch':'gateway-sync/'+('b'*40)+'/'+('c'*40),
                'state':'opened','sha':sha,'detailed_merge_status':'mergeable',
                'head_pipeline':{'id':8,'sha':sha,'project_id':19,
                                 'source':'merge_request_event','status':'success'}}
            obj.api=Mock(side_effect=lambda method,path: {'commit':{'id':'c'*40}}
                         if path.endswith('/repository/branches/main') else mr)
            obj.fetch=Mock(return_value=sha)
            obj.ancestor=Mock(return_value=True)
            obj.pure_sync=Mock(return_value='b'*40)
            obj.upstream_gate=Mock(return_value='clean')
            def pages(path):
                if path.endswith('/jobs'):
                    return [{'name':n,'status':'success'} for n in s.REQUIRED_JOBS]
                if path.endswith('/discussions'):return []
                raise AssertionError('Duplicate review must not be fetched')
            obj.pages=Mock(side_effect=pages)
            self.assertEqual(obj.handle_mr(mr,'b'*40,'c'*40)['status'],'would_merge')
            obj.upstream_gate.return_value='waiting_upstream_ci'
            self.assertEqual(obj.handle_mr(mr,'b'*40,'c'*40)['status'],'waiting_upstream_ci')
            obj.upstream_gate.return_value='clean'
            obj.pure_sync.return_value=None
            self.assertEqual(obj.handle_mr(mr,'b'*40,'c'*40)['status'],'blocked_non_sync_changes')

    def test_changed_target_does_not_modify_adopted_branch(self):
        with tempfile.TemporaryDirectory() as d:
            obj=s.Sync(Path(d))
            mr={'iid':4, 'source_project_id':19, 'target_project_id':19,
                'target_branch':'main', 'source_branch':s.ADOPT_BRANCH,
                'state':'opened', 'sha':'a'*40}
            obj.api=Mock(return_value=mr)
            obj.fetch=Mock(return_value='a'*40)
            obj.ancestor=Mock(return_value=False)
            obj.push_branch=Mock(side_effect=AssertionError('Never mutate adopted branch'))
            self.assertEqual(obj.handle_mr(mr,'b'*40,'c'*40)['status'],'blocked_target_changed')

class RealGit(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.root=Path(self.tmp.name)
        self.obj=s.Sync(self.root)
        subprocess.run(['git','init','--bare',str(self.obj.repo)],check=True,capture_output=True)

    def tearDown(self):self.tmp.cleanup()

    def commit(self, files, parents=()):
        rows=[]
        for name,text in sorted(files.items()):
            oid=self.obj.git('hash-object','-w','--stdin',stdin=text).stdout.strip()
            rows.append(f'100644 blob {oid}\t{name}\n')
        tree=self.obj.git('mktree',stdin=''.join(rows)).stdout.strip()
        args=[]
        for p in parents:args+=['-p',p]
        return self.obj.git('commit-tree',tree,*args,stdin='test\n').stdout.strip()

    def test_conflict_free_merge_preserves_both_histories_and_private_files(self):
        base=self.commit({'engine':'base'})
        internal=self.commit({'engine':'base','private-notes':'private'},[base])
        gateway=self.commit({'engine':'upstream'},[base])
        merged=self.obj.merged_commit(internal,gateway)
        self.assertTrue(self.obj.ancestor(internal,merged))
        self.assertTrue(self.obj.ancestor(gateway,merged))
        self.assertEqual(self.obj.git('show',merged+':private-notes').stdout,'private')
        self.assertEqual(self.obj.git('show',merged+':engine').stdout,'upstream')

    def test_conflicting_merge_stops_without_creating_a_ref(self):
        base=self.commit({'engine':'base\n'})
        internal=self.commit({'engine':'private\n'},[base])
        gateway=self.commit({'engine':'public\n'},[base])
        self.assertIsNone(self.obj.merged_commit(internal,gateway))
        self.assertEqual(self.obj.git('for-each-ref','--format=%(refname)').stdout,'')

    def test_pure_sync_rejects_extra_code_even_on_worker_branch(self):
        base=self.commit({'engine':'base'})
        internal=self.commit({'engine':'base','private-notes':'private'},[base])
        gateway=self.commit({'engine':'upstream'},[base])
        merged=self.obj.merged_commit(internal,gateway)
        mr={'iid':7,'source_project_id':19,'target_project_id':19,'target_branch':'main',
            'author':{'id':64},'source_branch':f'gateway-sync/{gateway}/{internal}'}
        self.assertEqual(self.obj.pure_sync(mr,merged,gateway,internal),gateway)
        edited=self.commit({'engine':'unreviewed edit','private-notes':'private'},[merged])
        self.assertIsNone(self.obj.pure_sync(mr,edited,gateway,internal))

    def test_status_is_durable_and_atomic(self):
        self.obj.status('blocked_review',mr=4)
        self.assertIn('blocked_review',(self.root/'status.json').read_text())
        self.assertFalse((self.root/'status.tmp').exists())

if __name__=='__main__':unittest.main()
