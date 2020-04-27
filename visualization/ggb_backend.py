#!/usr/bin/python3

import asyncio
from pyppeteer import launch

import base64
from PIL import Image
from io import BytesIO


BACKEND_LINK = 'file:/users/vadymka/research/geometry_rl/geogebra-math-apps-bundle-5-0-533-0/test.html'

class GeogebraPyppeteerBackend:
    def __init__(self, config):
        self.config = config

    async def evaluate_command_get_lables(self, js_cmd):
        cmd = '''window.ggbApplet.evalCommandGetLabels("''' + js_cmd +'''")'''
        res = await self.page.evaluate(cmd)
        return res


    async def check_condition(self, cond):
        # print('bool_value = ' + cond)
        await self.evaluate_command('bool_value = ' + cond)
        ans = await self.get_value('bool_value')
        return ans


    async def set_undo_point(self):
        cmd = '''window.ggbApplet.setUndoPoint()'''
        await self.page.evaluate(cmd)


    async def undo(self):
        cmd = '''window.ggbApplet.undo()'''
        await self.page.evaluate(cmd)


    async def is_defined(self, obj):
        cmd = '''window.ggbApplet.isDefined("''' + obj + '''")'''
        return await self.page.evaluate(cmd)


    async def delete(self, obj):
        # print('delete:', obj)
        cmd = '''window.ggbApplet.deleteObject("''' + obj + '''")'''
        await self.page.evaluate(cmd)


    async def evaluate_command(self, js_cmd):
        cmd = '''window.ggbApplet.evalCommand("''' + js_cmd +'''")'''
        res = await self.page.evaluate(cmd)
        return res


    async def get_value(self, obj):
        cmd = '''window.ggbApplet.getValue("''' + obj +'''")'''
        # print(cmd)
        res = await self.page.evaluate(cmd)
        return res


    async def get_type(self, obj):
        cmd = '''window.ggbApplet.getObjectType("''' + obj +'''")'''
        res = await self.page.evaluate(cmd)
        return res


    async def get_objects(self, obj_type=''):
        cmd = '''window.ggbApplet.getAllObjectNames("''' + obj_type +'''")'''
        res = await self.page.evaluate(cmd)
        return res


    async def are_equal(self, o1, o2):
        cmd = 'AreEqual({0},{1})'.format(o1, o2)
        # print('get_value:', cmd)
        value = await self.get_value(cmd)
        return value > 0.5


    async def run(self):
        self.browser = await launch(
            headless=self.config['headless'],
            args=['--window-size=3000,2000','--no-sandbox']
        )
        self.page = await self.browser.newPage()
        await self.page.setViewport({"width": 1200, "height": 1000, "deviceScaleFactor": 2})
        await self.page.goto(BACKEND_LINK)
        await self.page.waitForFunction('window.ggbApplet != null');
        self.applet = await self.page.evaluate('window.ggbApplet')
        #???
        print('ggb backend started')


    async def stop(self):
        await self.browser.close()
        print('ggb backend closed')


    async def applet_screenshot(self, path='./data/example.png'):
        print('save screenshot to:', path)
        cmd = '''window.ggbApplet.getPNGBase64(1, false, 720)'''

        data = await self.page.evaluate(cmd)
        im = Image.open(BytesIO(base64.b64decode(data)))
        im.save(path, 'PNG')


    async def screenshot(self, path='./data/example.png'):
        print('save screenshot to:', path)
        await self.page.screenshot({'path': path, 'quality':100})


    # def async def load_configuration(self):
    async def init(self):
        while False:
            try:
                print(await e.ggb.are_equal('A', 'A'))
                break
            except:
                continue


        await self.evaluate_command('A = Point({1.1, 5.45})');
        await self.evaluate_command('B = Point({5, 0})');
        await self.evaluate_command('C = Point({-2, -2})');

        await self.evaluate_command('la = Line(B, C)');
        await self.evaluate_command('lb = Line(A, C)');
        await self.evaluate_command('lc = Line(A, B)');

        await self.evaluate_command('ha = PerpendicularLine(A, la)');
        await self.evaluate_command('hb = PerpendicularLine(B, lb)');
        await self.evaluate_command('hc = PerpendicularLine(C, lc)');

        await self.evaluate_command('H = Intersect(ha, hb)');
        await self.evaluate_command('Ha = Intersect(la, ha)');
        await self.evaluate_command('Hb = Intersect(lb, hb)');
        await self.evaluate_command('Hc = Intersect(lc, hc)');

        await self.evaluate_command('W = Circle(A, B, C)');
        await self.evaluate_command('O = Center(W)');
        await self.evaluate_command('e = Line(H, O)');

        await self.evaluate_command('Ma = Midpoint(B, C)');
        await self.evaluate_command('Mb = Midpoint(A, C)');
        await self.evaluate_command('Mc = Midpoint(A, B)');

        await self.evaluate_command('l1 = Line(H, Ma)');
        ans = await self.evaluate_command_get_lables('X = Intersect(W, ha)');
        print(ans)


async def test_backend():
    gb = GeogebraPyppeteerBackend({'headless' : True})

    await gb.run()
    await gb.init()

    print(await gb.get_objects('point'))
    print(await gb.get_objects('line'))
    print(await gb.get_objects('circle'))

    # print(await gb.evaluate_command('l '))

    await gb.screenshot()
    await gb.stop()

    # display(Image('./data/example.png'))

def main():
    loop = asyncio.get_event_loop()
    loop.run_until_complete(test_backend())
    loop.close()

if __name__ == "__main__":
    main()


