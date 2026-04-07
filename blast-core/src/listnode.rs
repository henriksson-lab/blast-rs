//! Rust equivalent of the ListNode linked list from ncbi_std.c.
//! In Rust, we use Vec instead of a linked list for most purposes.

/// A simple linked list node equivalent.
/// In the Rust rewrite, most uses of ListNode will be replaced with Vec.
#[derive(Debug)]
pub struct ListNode<T> {
    pub choice: u8,
    pub data: T,
}

/// A list of tagged data items (replaces the C ListNode linked list).
pub type NodeList<T> = Vec<ListNode<T>>;

/// Create a new empty list.
pub fn new_list<T>() -> NodeList<T> {
    Vec::new()
}

/// Add an item to the end of the list.
pub fn add_pointer<T>(list: &mut NodeList<T>, choice: u8, data: T) {
    list.push(ListNode { choice, data });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_and_access() {
        let mut list: NodeList<String> = new_list();
        add_pointer(&mut list, 1, "hello".to_string());
        add_pointer(&mut list, 2, "world".to_string());
        assert_eq!(list.len(), 2);
        assert_eq!(list[0].choice, 1);
        assert_eq!(list[0].data, "hello");
        assert_eq!(list[1].choice, 2);
    }
}
