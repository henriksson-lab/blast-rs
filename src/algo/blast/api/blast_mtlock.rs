use std::any::Any;
use std::fmt;
use std::sync::{Arc, Mutex};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EMTLock {
    Lock,
    LockRead,
    Unlock,
    TryLock,
    TryLockRead,
}

pub type FMTLockHandler = fn(&mut dyn Any, EMTLock) -> i32;
pub type FMTLockCleanup = fn(Box<dyn Any + Send>);

pub struct MTLock {
    inner: Arc<MTLockInner>,
}

struct MTLockInner {
    data: Mutex<Option<Box<dyn Any + Send>>>,
    handler: Option<FMTLockHandler>,
    cleanup: Option<FMTLockCleanup>,
}

impl Clone for MTLock {
    fn clone(&self) -> Self {
        Self {
            inner: Arc::clone(&self.inner),
        }
    }
}

impl fmt::Debug for MTLock {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MTLock")
            .field("strong_count", &Arc::strong_count(&self.inner))
            .field("has_handler", &self.inner.handler.is_some())
            .field("has_cleanup", &self.inner.cleanup.is_some())
            .finish()
    }
}

impl Drop for MTLockInner {
    fn drop(&mut self) {
        if let Some(cleanup) = self.cleanup {
            if let Some(data) = self.data.lock().unwrap().take() {
                cleanup(data);
            }
        }
    }
}

// --- function: MT_LOCK_Create -> mt_lock_create ---
pub fn mt_lock_create(
    data: Option<Box<dyn Any + Send>>,
    handler: Option<FMTLockHandler>,
    cleanup: Option<FMTLockCleanup>,
) -> Option<MTLock> {
    Some(MTLock {
        inner: Arc::new(MTLockInner {
            data: Mutex::new(data),
            handler,
            cleanup,
        }),
    })
}

// --- function: MT_LOCK_AddRef -> mt_lock_add_ref ---
pub fn mt_lock_add_ref(lk: Option<&MTLock>) -> Option<MTLock> {
    let lk = lk?;
    let clone = lk.clone();
    if lk.inner.handler.is_some() {
        let _ = mt_lock_do(Some(lk), EMTLock::Lock);
        let _ = mt_lock_do(Some(lk), EMTLock::Unlock);
    }
    Some(clone)
}

// --- function: MT_LOCK_Delete -> mt_lock_delete ---
pub fn mt_lock_delete(lk: Option<MTLock>) -> Option<MTLock> {
    drop(lk);
    None
}

// --- function: MT_LOCK_DoInternal / MT_LOCK_Do -> mt_lock_do ---
pub fn mt_lock_do(lk: Option<&MTLock>, how: EMTLock) -> i32 {
    let Some(lk) = lk else {
        return -1;
    };
    let Some(handler) = lk.inner.handler else {
        return -1;
    };
    let mut data = lk.inner.data.lock().unwrap();
    if let Some(data) = data.as_deref_mut() {
        handler(data, how)
    } else {
        let mut data = ();
        handler(&mut data, how)
    }
}
