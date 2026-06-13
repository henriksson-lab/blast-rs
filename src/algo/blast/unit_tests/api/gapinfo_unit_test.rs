// Upstream source:
// - ncbi-blast-2.17.0+-src/c++/src/algo/blast/unit_tests/api/gapinfo_unit_test.cpp

#[cfg(test)]
mod tests {
    use crate::algo::blast::core::gapinfo::{
        gap_edit_script_delete, gap_edit_script_dup, gap_edit_script_new, GapAlignOpType,
    };

    #[test]
    fn test_gap_edit_script_new() {
        let size = 3;

        let esp = gap_edit_script_new(size).expect("positive edit script size should allocate");
        assert_eq!(size, esp.size);
        assert!(gap_edit_script_delete(Some(esp)).is_none());

        assert!(gap_edit_script_new(-1).is_none());
    }

    #[test]
    fn test_gap_edit_script_dup() {
        let size = 3;
        let nums = [7, 11, 13];
        let op_types = [
            GapAlignOpType::Sub,
            GapAlignOpType::Del,
            GapAlignOpType::Ins,
        ];
        let mut esp =
            gap_edit_script_new(size).expect("positive edit script size should allocate");

        for index in 0..size as usize {
            esp.num[index] = nums[index];
            esp.op_type[index] = op_types[index];
        }

        let esp_dup = gap_edit_script_dup(Some(&esp)).expect("valid edit script should duplicate");
        assert_eq!(size, esp_dup.size);

        for index in 0..size as usize {
            assert_eq!(nums[index], esp_dup.num[index]);
            assert_eq!(op_types[index], esp_dup.op_type[index]);
        }

        assert!(gap_edit_script_delete(Some(esp)).is_none());
        assert!(gap_edit_script_delete(Some(esp_dup)).is_none());
    }
}
