import json
import socket
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import threading
import webbrowser

from deconvolution import analyze

ROOT = Path(__file__).resolve().parent
SAMPLES = {path.name: path for path in (ROOT / 'samples').glob('*.txt')}


def port_in_use(port):
    """Windows 上 SO_REUSEADDR 允许两个进程绑同一端口，先探测再绑定。"""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.settimeout(0.3)
        return probe.connect_ex(('127.0.0.1', port)) == 0


class Handler(BaseHTTPRequestHandler):
    def send_json(self, data, status=200):
        encoded = json.dumps(data, ensure_ascii=False).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Content-Length', str(len(encoded)))
        self.end_headers()
        self.wfile.write(encoded)

    def do_GET(self):
        if self.path == '/api/samples':
            self.send_json(list(SAMPLES))
            return
        if self.path in ('/', '/index.html'):
            encoded = (ROOT / 'index.html').read_bytes()
            self.send_response(200)
            self.send_header('Content-Type', 'text/html; charset=utf-8')
            self.send_header('Content-Length', str(len(encoded)))
            self.end_headers()
            self.wfile.write(encoded)
            return
        self.send_json({'error': 'Not found'}, 404)

    def do_POST(self):
        if self.path != '/api/analyze':
            self.send_json({'error': 'Not found'}, 404)
            return
        if self.headers.get('Origin') not in (None, f'http://127.0.0.1:{self.server.server_port}', f'http://localhost:{self.server.server_port}'):
            self.send_json({'error': '不接受其他网站发来的分析请求。'}, 403)
            return
        try:
            length = int(self.headers.get('Content-Length', 0))
            if not 0 < length <= 25_000_000:
                raise ValueError('请选择 25 MB 以内的两列文本文件。')
            try:
                request = json.loads(self.rfile.read(length))
            except json.JSONDecodeError as error:
                raise ValueError('请求体不是合法 JSON。') from error
            if not isinstance(request, dict):
                raise ValueError('请求体必须是 JSON 对象。')
            sample = request.get('sample')
            if sample:
                if not isinstance(sample, str) or sample not in SAMPLES:
                    raise ValueError('示例文件不存在。')
                text = SAMPLES[sample].read_text(encoding='utf-8-sig')
                filename = sample
            else:
                text = request.get('text', '')
                if not isinstance(text, str):
                    raise ValueError('上传内容必须是两列文本。')
                filename = str(request.get('filename', 'spectrum.txt'))
            settings = request.get('settings', {})
            if not isinstance(settings, dict):
                raise ValueError('settings 必须是 JSON 对象。')
            result = analyze(text, filename=filename, **settings)
            self.send_json(result)
        except (ValueError, TypeError, KeyError) as error:
            self.send_json({'error': str(error)}, 400)
        except Exception as error:
            self.send_json({'error': f'分析未完成：{type(error).__name__}: {error}'}, 500)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', type=int, default=8765)
    parser.add_argument('--no-browser', action='store_true')
    args = parser.parse_args()
    if port_in_use(args.port):
        raise SystemExit(f'端口 {args.port} 已被占用。可改用 python server.py --port {args.port + 1}。')
    try:
        server = ThreadingHTTPServer(('127.0.0.1', args.port), Handler)
    except OSError as error:
        raise SystemExit(f'端口 {args.port} 不可用。可改用 python server.py --port {args.port + 1}。\n{error}')
    print(f'HisMS 已启动：http://127.0.0.1:{args.port}  关闭窗口或按 Ctrl+C 停止。', flush=True)
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(f'http://127.0.0.1:{args.port}')).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
