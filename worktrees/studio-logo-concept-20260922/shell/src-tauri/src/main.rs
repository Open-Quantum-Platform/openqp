#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

// The desktop UI is bundled into Tauri. Its Python sidecar receives requests
// over stdio, so the installed app never needs a loopback HTTP server.

use std::{
    collections::HashMap,
    sync::{
        atomic::{AtomicU64, Ordering},
        mpsc, Mutex,
    },
    time::Duration,
};

use serde::{Deserialize, Serialize};
use tauri::Manager;
use tauri_plugin_shell::{
    process::{CommandChild, CommandEvent},
    ShellExt,
};

#[derive(Deserialize, Serialize)]
struct RpcRequest {
    method: String,
    path: String,
    headers: HashMap<String, String>,
    body: String,
}

#[derive(Clone, Deserialize, Serialize)]
struct RpcResponse {
    status: u16,
    headers: HashMap<String, String>,
    body: String,
}

#[derive(Deserialize)]
struct RpcMessage {
    id: u64,
    result: Option<RpcResponse>,
    error: Option<String>,
}

struct Backend {
    child: Mutex<Option<CommandChild>>,
    pending: Mutex<HashMap<u64, mpsc::Sender<Result<RpcResponse, String>>>>,
    next_id: AtomicU64,
}

impl Backend {
    fn call(&self, request: RpcRequest) -> Result<RpcResponse, String> {
        let id = self.next_id.fetch_add(1, Ordering::Relaxed);
        let (sender, receiver) = mpsc::channel();
        self.pending
            .lock()
            .map_err(|_| "backend request state is unavailable".to_string())?
            .insert(id, sender);

        let payload = serde_json::json!({ "id": id, "request": request }).to_string() + "\n";
        let write_result = self
            .child
            .lock()
            .map_err(|_| "backend process state is unavailable".to_string())?
            .as_mut()
            .ok_or_else(|| "backend process is not running".to_string())?
            .write(payload.as_bytes());
        if let Err(error) = write_result {
            self.pending
                .lock()
                .ok()
                .and_then(|mut pending| pending.remove(&id));
            return Err(error.to_string());
        }

        match receiver.recv_timeout(Duration::from_secs(180)) {
            Ok(result) => result,
            Err(_) => {
                self.pending
                    .lock()
                    .ok()
                    .and_then(|mut pending| pending.remove(&id));
                Err("backend request timed out".to_string())
            }
        }
    }

    fn resolve(&self, message: RpcMessage) {
        let result = match (message.result, message.error) {
            (Some(response), _) => Ok(response),
            (_, Some(error)) => Err(error),
            _ => Err("backend returned an invalid response".to_string()),
        };
        if let Ok(mut pending) = self.pending.lock() {
            if let Some(sender) = pending.remove(&message.id) {
                let _ = sender.send(result);
            }
        }
    }
}

#[tauri::command]
fn backend_call(
    state: tauri::State<'_, Backend>,
    request: RpcRequest,
) -> Result<RpcResponse, String> {
    state.call(request)
}

fn sidecar_environment(app: &tauri::App) -> HashMap<String, String> {
    let mut environment = HashMap::new();
    if let Ok(resources) = app.path().resource_dir() {
        environment.insert(
            "OQP_STUDIO_RESOURCES".to_string(),
            resources.to_string_lossy().to_string(),
        );
    }
    environment
}

fn main() {
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .manage(Backend {
            child: Mutex::new(None),
            pending: Mutex::new(HashMap::new()),
            next_id: AtomicU64::new(1),
        })
        .setup(|app| {
            let (receiver, child) = app
                .shell()
                .sidecar("oqp-studio-backend")?
                .arg("--stdio")
                .envs(sidecar_environment(app))
                .spawn()?;
            *app.state::<Backend>().child.lock().unwrap() = Some(child);

            let handle = app.handle().clone();
            std::thread::spawn(move || {
                let mut receiver = receiver;
                while let Some(event) = tauri::async_runtime::block_on(receiver.recv()) {
                    if let CommandEvent::Stdout(bytes) = event {
                        let Ok(line) = String::from_utf8(bytes) else {
                            continue;
                        };
                        if let Ok(message) = serde_json::from_str::<RpcMessage>(line.trim()) {
                            handle.state::<Backend>().resolve(message);
                        }
                    }
                }
            });
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![backend_call])
        .build(tauri::generate_context!())
        .expect("error while building OQP Studio")
        .run(|app, event| {
            if let tauri::RunEvent::Exit = event {
                if let Some(child) = app.state::<Backend>().child.lock().unwrap().take() {
                    let _ = child.kill();
                }
            }
        });
}
