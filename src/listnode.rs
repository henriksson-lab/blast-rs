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

/// Rust equivalent of `ListNodeAdd` (`ncbi_std.c:77`).
pub fn list_node_add<T: Default>(list: &mut NodeList<T>) -> &mut ListNode<T> {
    list.push(ListNode {
        choice: 0,
        data: T::default(),
    });
    list.last_mut().expect("just pushed a ListNode")
}

/// Rust equivalent of `ListNodeAddPointer` (`ncbi_std.c:99`).
pub fn list_node_add_pointer<T>(list: &mut NodeList<T>, choice: u8, data: T) -> &mut ListNode<T> {
    list.push(ListNode { choice, data });
    list.last_mut().expect("just pushed a ListNode")
}

/// Rust equivalent of `ListNodeCopyStr` (`ncbi_std.c:124`).
pub fn list_node_copy_str<'a>(
    list: &'a mut NodeList<String>,
    choice: u8,
    str_: Option<&str>,
) -> Option<&'a mut ListNode<String>> {
    let str_ = str_?;
    Some(list_node_add_pointer(list, choice, str_.to_owned()))
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

    #[test]
    fn test_c_shape_list_node_helpers() {
        let mut list: NodeList<String> = new_list();
        let node = list_node_add(&mut list);
        node.choice = 7;
        node.data = "defaulted".to_owned();
        assert_eq!(list[0].choice, 7);
        assert_eq!(list[0].data, "defaulted");

        let node = list_node_add_pointer(&mut list, 3, "owned".to_owned());
        assert_eq!(node.choice, 3);
        assert_eq!(node.data, "owned");

        let node = list_node_copy_str(&mut list, 4, Some("copied")).unwrap();
        assert_eq!(node.choice, 4);
        assert_eq!(node.data, "copied");
        assert!(list_node_copy_str(&mut list, 9, None).is_none());
        assert_eq!(list.len(), 3);
    }
}
