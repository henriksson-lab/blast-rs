use crate::algo::blast::api::local_db_adapter::CLocalDbAdapter;
use crate::algo::blast::api::psiblast_impl::CPsiBlastImpl;

pub struct CPsiBlast {
    pub subject: Option<Box<CLocalDbAdapter>>,
    pub impl_: *mut CPsiBlastImpl,
}
