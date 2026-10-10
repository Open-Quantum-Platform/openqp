const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const assert=require('node:assert/strict');
(async()=>{
 const browser=await chromium.launch({headless:true,executablePath:process.env.CHROME_EXECUTABLE,args:['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 try {
 const page=await browser.newPage({viewport:{width:1400,height:950}});
 const session=await page.context().newCDPSession(page); await session.send('Emulation.setCPUThrottlingRate',{rate:6});
 const errors=[];page.on('pageerror',e=>errors.push(e.message));
 await page.route('**/api/health',async route=>{await new Promise(r=>setTimeout(r,6000));await route.continue();});
 await page.goto(process.env.STUDIO_TEST_URL || 'http://127.0.0.1:15173');
 await page.waitForFunction(()=>document.getElementById('health').textContent.includes('1s'));
 assert.equal(await page.locator('#health').getAttribute('aria-busy'),'true');
 await page.screenshot({path:'startup-clock.png'});
 await page.locator('[data-tab="method"]').click();
 assert.equal(await page.locator('#panel-method').getAttribute('aria-hidden'),'false');
 await page.waitForFunction(()=>document.getElementById('health').textContent.includes('ready'));
 await page.locator('[data-tab="art"]').click();
 const frame=page.frameLocator('#artFrame');
 await frame.locator('#progress').filter({hasText:'Interactive'}).waitFor();
 assert.equal(await frame.locator('#renderMode').inputValue(),'preview');
 assert.equal(await frame.locator('#renderMode option[value="ray"]').evaluate(e=>e.disabled),true);
 console.log('gpu', await frame.locator('canvas').evaluate(c=>{const g=c.getContext('webgl2'), e=g.getExtension('WEBGL_debug_renderer_info');return {name:e?g.getParameter(e.UNMASKED_RENDERER_WEBGL):'unknown',userAgent:navigator.userAgent,cores:navigator.hardwareConcurrency};}));
 await page.screenshot({path:'art-interactive.png'});
 const dims=await frame.locator('canvas').evaluate(c=>({w:c.width,h:c.height}));
 assert.ok(dims.w*dims.h<=1000001);
 await frame.locator('body').evaluate(()=>{
   window.artMessages=[];
   window.addEventListener('message',e=>{
     if(e.data?.type==='oqp-art-scene'||e.data?.type==='oqp-art-active')window.artMessages.push(e.data.type);
   });
 });
 await page.locator('[data-tab="method"]').click();
 await page.waitForTimeout(250);
 assert.equal(await frame.locator('#renderMode').inputValue(),'preview');
 await page.locator('[data-tab="art"]').click();
 await frame.locator('#progress').filter({hasText:'Interactive'}).waitFor();
 assert.deepEqual(await frame.locator('body').evaluate(()=>window.artMessages.slice(-2)),['oqp-art-scene','oqp-art-active']);
 // A rejected large structure must never redisplay the previous orbital.
 const values=Array.from({length:27},(_,i)=>Math.floor(i/9)-1).join(' ');
 const cube='test\norbital\n0 0 0 0\n3 1 0 0\n3 0 1 0\n3 0 0 1\n'+values;
 await page.evaluate(cube=>document.querySelector('#artFrame').contentWindow.postMessage(
   {type:'oqp-art-scene',scene:{atoms:[],cube:'data:text/plain,'+encodeURIComponent(cube),iso:0.1}},location.origin),cube);
 await page.waitForTimeout(1200);
 assert.match(await frame.locator('#progress').textContent(),/Interactive preview/);
 // While a newer cube downloads, neither tracing nor exporting stale art is allowed.
 let releaseCube;
 const cubeGate=new Promise(resolve=>{releaseCube=resolve;});
 await page.route('**/delayed-art.cube',async route=>{
   await cubeGate;
   await route.fulfill({status:200,contentType:'text/plain',body:cube});
 });
 await page.evaluate(()=>document.querySelector('#artFrame').contentWindow.postMessage(
   {type:'oqp-art-scene',scene:{atoms:[],cube:location.origin+'/delayed-art.cube',iso:0.1}},location.origin));
 await frame.locator('#progress').filter({hasText:'Loading artwork'}).waitFor();
 assert.equal(await frame.locator('#renderMode').isDisabled(),true);
 assert.equal(await frame.locator('#save').isDisabled(),true);
 // Even a dispatched change must be refused while loading (also safe on a real GPU).
 await frame.locator('#renderMode').evaluate(e=>{e.value='ray';e.dispatchEvent(new Event('change'));});
 assert.equal(await frame.locator('#renderMode').inputValue(),'preview');
 releaseCube();
 await page.waitForFunction(()=>!document.querySelector('#artFrame').contentDocument.querySelector('#renderMode').disabled);
 await frame.locator('#progress').filter({hasText:'Interactive preview'}).waitFor();
 assert.equal(await frame.locator('#save').isDisabled(),false);
 await frame.locator('#surfacePhases').selectOption('positive');
 await page.waitForTimeout(250);
 await page.evaluate(()=>document.querySelector('#artFrame').contentWindow.postMessage(
   {type:'oqp-art-scene',scene:{atoms:Array.from({length:301},()=>['C',0,0,0])}},location.origin));
 await frame.locator('#progress').filter({hasText:'up to 300 atoms'}).waitFor();
 await frame.locator('#surfacePhases').selectOption('negative');
 await page.waitForTimeout(250);
 assert.match(await frame.locator('#progress').textContent(),/up to 300 atoms/);
 assert.deepEqual(errors,[]);
 console.log(JSON.stringify({checks:['startup clock increments','UI navigation while backend waits','ready state','software GPU safe mode','pixel budget','hidden Art stops tracing','Art reopens with one deferred scene activation','cube loading prevents stale tracing/export','oversized rejection remains clear after material changes'],canvas:dims,pageErrors:errors}));
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exit(1)});
