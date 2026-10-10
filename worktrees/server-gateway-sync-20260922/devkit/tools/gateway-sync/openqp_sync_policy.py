"""Recognize only the dedicated inbound worker's MRs and the approved initial MR."""
import re

INITIAL_GATEWAY = '255aef4ac87cae7c18f70d987ebad9b052f76481'
INITIAL_SHA = '2dd8b09efc2c013667256bd06bf67305d4f6ad19'
INITIAL_BRANCH = 'codex/gateway-internal-sync-20260922'
INITIAL_DOC = 'devkit/docs/gateway-inbound-sync.md'
INITIAL_DOC_BLOB = 'a1b1d89343b20f5f355a3585b663bf541495d191'


def sync_gateway(project_id, mr):
    if (project_id != 19 or mr.get('source_project_id') != 19
            or mr.get('target_project_id') != 19 or mr.get('target_branch') != 'main'):
        return None
    author = mr.get('author', {}).get('id')
    branch = mr.get('source_branch', '')
    if (mr.get('iid') == 4 and author == 61 and branch == INITIAL_BRANCH
            and mr.get('sha') == INITIAL_SHA):
        return INITIAL_GATEWAY
    match = re.fullmatch(r'gateway-sync/([0-9a-f]{40})/([0-9a-f]{40})', branch)
    if author == 64 and match:
        return match[1]
    return None
