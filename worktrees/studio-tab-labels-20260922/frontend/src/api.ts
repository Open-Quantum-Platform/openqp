type RpcResponse = {
  status: number;
  headers: Record<string, string>;
  body: string;
};

type TauriWindow = Window & {
  __TAURI__?: { core?: { invoke: (command: string, args: unknown) => Promise<RpcResponse> } };
};

function base64(bytes: Uint8Array): string {
  let text = "";
  for (const byte of bytes) text += String.fromCharCode(byte);
  return btoa(text);
}

function fromBase64(value: string): ArrayBuffer {
  const text = atob(value);
  const bytes = new Uint8Array(text.length);
  for (let index = 0; index < text.length; index += 1) bytes[index] = text.charCodeAt(index);
  return bytes.buffer;
}

function apiPath(input: RequestInfo | URL): string | null {
  const url = new URL(input instanceof Request ? input.url : input, window.location.origin);
  return url.pathname.startsWith("/api/") ? `${url.pathname}${url.search}` : null;
}

/** Stop waiting when a caller cancels; the sidecar may still finish its work. */
export function withAbort<T>(operation: Promise<T>, signal: AbortSignal): Promise<T> {
  return new Promise<T>((resolve, reject) => {
    const abort = () => reject(signal.reason ?? new DOMException("Aborted", "AbortError"));
    if (signal.aborted) {
      operation.catch(() => {});
      abort();
      return;
    }
    signal.addEventListener("abort", abort, { once: true });
    operation.then(resolve, reject).finally(() => signal.removeEventListener("abort", abort));
  });
}

/** Route API calls through Tauri IPC in the desktop bundle, while keeping Vite
 * development on ordinary fetch so the browser loop remains immediate. */
export function installApiFetch(): void {
  const invoke = (window as TauriWindow).__TAURI__?.core?.invoke;
  if (!invoke) return;

  const webFetch = window.fetch.bind(window);
  window.fetch = async (input: RequestInfo | URL, init?: RequestInit): Promise<Response> => {
    const path = apiPath(input);
    if (!path) return webFetch(input, init);

    const request = new Request(input, init);
    const body = request.method === "GET" || request.method === "HEAD"
      ? ""
      : base64(new Uint8Array(await request.arrayBuffer()));
    request.signal.throwIfAborted();
    const result = await withAbort(invoke("backend_call", {
      request: {
        method: request.method,
        path,
        headers: Object.fromEntries(request.headers.entries()),
        body,
      },
    }), request.signal);
    return new Response(fromBase64(result.body), {
      status: result.status,
      headers: result.headers,
    });
  };
}
