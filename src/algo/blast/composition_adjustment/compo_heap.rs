use std::ffi::c_void;

const COMPO_INTENSE_DEBUG: bool = false;
const HEAP_INITIAL_CAPACITY: i32 = 100;
const HEAP_RESIZE_FACTOR: f64 = 1.5;
const HEAP_MIN_RESIZE: i32 = 100;

#[derive(Clone, Copy, Debug)]
pub struct BlastCompoHeapRecord {
    pub best_evalue: f64,
    pub best_score: i32,
    pub subject_index: i32,
    pub these_alignments: *mut c_void,
}

#[derive(Clone, Debug)]
pub struct BlastCompoHeap {
    pub n: i32,
    pub capacity: i32,
    pub heap_threshold: i32,
    pub ecutoff: f64,
    pub worst_evalue: f64,
    pub array: Option<Vec<BlastCompoHeapRecord>>,
    pub heap_array: Option<Vec<BlastCompoHeapRecord>>,
}

fn cmp<T: PartialOrd>(a: T, b: T) -> i32 {
    if a > b {
        1
    } else if a < b {
        -1
    } else {
        0
    }
}

fn s_compo_heap_record_compare(
    place1: &BlastCompoHeapRecord,
    place2: &BlastCompoHeapRecord,
) -> bool {
    let mut result = cmp(place1.best_evalue, place2.best_evalue);
    if result == 0 {
        result = cmp(place2.best_score, place1.best_score);
        if result == 0 {
            result = cmp(place2.subject_index, place1.subject_index);
        }
    }
    result > 0
}

fn s_compo_heap_record_swap(
    record1: &mut BlastCompoHeapRecord,
    record2: &mut BlastCompoHeapRecord,
) {
    std::mem::swap(record1, record2);
}

fn s_compo_heap_is_valid(heap_array: &[BlastCompoHeapRecord], i: i32, n: i32) -> bool {
    let left = 2 * i;
    let right = 2 * i + 1;

    if right <= n {
        return !s_compo_heap_record_compare(&heap_array[right as usize], &heap_array[i as usize])
            && s_compo_heap_is_valid(heap_array, right, n);
    }
    if left <= n {
        return !s_compo_heap_record_compare(&heap_array[left as usize], &heap_array[i as usize])
            && s_compo_heap_is_valid(heap_array, left, n);
    }
    true
}

fn s_compo_heapify_down(heap_array: &mut [BlastCompoHeapRecord], top: i32, n: i32) {
    let mut largest = top;
    loop {
        let i = largest;
        let left = 2 * i;
        let right = 2 * i + 1;

        if left <= n
            && s_compo_heap_record_compare(&heap_array[left as usize], &heap_array[i as usize])
        {
            largest = left;
        } else {
            largest = i;
        }

        if right <= n
            && s_compo_heap_record_compare(
                &heap_array[right as usize],
                &heap_array[largest as usize],
            )
        {
            largest = right;
        }

        if largest != i {
            let a = i as usize;
            let b = largest as usize;
            if a < b {
                let (left_part, right_part) = heap_array.split_at_mut(b);
                s_compo_heap_record_swap(&mut left_part[a], &mut right_part[0]);
            } else {
                let (left_part, right_part) = heap_array.split_at_mut(a);
                s_compo_heap_record_swap(&mut right_part[0], &mut left_part[b]);
            }
        } else {
            break;
        }
    }

    if COMPO_INTENSE_DEBUG {
        assert!(s_compo_heap_is_valid(heap_array, top, n));
    }
}

fn s_compo_heapify_up(heap_array: &mut [BlastCompoHeapRecord], mut i: i32) {
    let mut parent = i / 2;
    while parent >= 1
        && s_compo_heap_record_compare(&heap_array[i as usize], &heap_array[parent as usize])
    {
        let child = i as usize;
        let parent_idx = parent as usize;
        let (left, right) = heap_array.split_at_mut(child);
        s_compo_heap_record_swap(&mut right[0], &mut left[parent_idx]);

        i = parent;
        parent /= 2;
    }

    if COMPO_INTENSE_DEBUG {
        assert!(s_compo_heap_is_valid(heap_array, 1, i));
    }
}

fn s_convert_to_heap(self_: &mut BlastCompoHeap) {
    if let Some(array) = self_.array.take() {
        self_.heap_array = Some(array);
        let n = self_.n;
        if let Some(heap_array) = self_.heap_array.as_mut() {
            for i in (1..=(n / 2)).rev() {
                s_compo_heapify_down(heap_array, i, n);
            }
            if COMPO_INTENSE_DEBUG {
                assert!(s_compo_heap_is_valid(heap_array, 1, self_.n));
            }
        }
    }
}

pub fn blast_compo_heap_would_insert(
    self_: &mut BlastCompoHeap,
    evalue: f64,
    score: i32,
    subject_index: i32,
) -> bool {
    if self_.n < self_.heap_threshold || evalue <= self_.ecutoff || evalue < self_.worst_evalue {
        true
    } else {
        if self_.heap_array.is_none() {
            s_convert_to_heap(self_);
        }

        let heap_record = BlastCompoHeapRecord {
            best_evalue: evalue,
            best_score: score,
            subject_index,
            these_alignments: std::ptr::null_mut(),
        };

        let heap_array = self_
            .heap_array
            .as_ref()
            .expect("heap array is initialized");
        s_compo_heap_record_compare(&heap_array[1], &heap_record)
    }
}

fn s_comp_heap_record_insert_at_end(
    array: &mut Vec<BlastCompoHeapRecord>,
    length: &mut i32,
    capacity: &mut i32,
    alignments: *mut c_void,
    evalue: f64,
    score: i32,
    subject_index: i32,
) -> i32 {
    if *length >= *capacity {
        let new_capacity =
            (HEAP_MIN_RESIZE + *capacity).max((HEAP_RESIZE_FACTOR * (*capacity as f64)) as i32);
        array.resize(
            (new_capacity + 1) as usize,
            BlastCompoHeapRecord {
                best_evalue: 0.0,
                best_score: 0,
                subject_index: 0,
                these_alignments: std::ptr::null_mut(),
            },
        );
        *capacity = new_capacity;
    }

    *length += 1;
    array[*length as usize] = BlastCompoHeapRecord {
        best_evalue: evalue,
        best_score: score,
        subject_index,
        these_alignments: alignments,
    };

    0
}

pub fn blast_compo_heap_insert(
    self_: &mut BlastCompoHeap,
    alignments: *mut c_void,
    evalue: f64,
    score: i32,
    subject_index: i32,
    discarded_alignments: &mut *mut c_void,
) -> i32 {
    *discarded_alignments = std::ptr::null_mut();

    if self_.array.is_some() && self_.n >= self_.heap_threshold {
        s_convert_to_heap(self_);
    }

    if let Some(array) = self_.array.as_mut() {
        let status = s_comp_heap_record_insert_at_end(
            array,
            &mut self_.n,
            &mut self_.capacity,
            alignments,
            evalue,
            score,
            subject_index,
        );
        if status != 0 {
            return -1;
        }
        if self_.worst_evalue < evalue {
            self_.worst_evalue = evalue;
        }
    } else {
        let heap_array = self_
            .heap_array
            .as_mut()
            .expect("heap array is initialized");
        if self_.n < self_.heap_threshold
            || (evalue <= self_.ecutoff && self_.worst_evalue <= self_.ecutoff)
        {
            let status = s_comp_heap_record_insert_at_end(
                heap_array,
                &mut self_.n,
                &mut self_.capacity,
                alignments,
                evalue,
                score,
                subject_index,
            );
            if status != 0 {
                return -1;
            }
            s_compo_heapify_up(heap_array, self_.n);
        } else {
            let heap_record = BlastCompoHeapRecord {
                best_evalue: evalue,
                best_score: score,
                subject_index,
                these_alignments: alignments,
            };

            if s_compo_heap_record_compare(&heap_array[1], &heap_record) {
                *discarded_alignments = heap_array[1].these_alignments;
                heap_array[1] = heap_record;
            } else {
                *discarded_alignments = heap_record.these_alignments;
            }
            s_compo_heapify_down(heap_array, 1, self_.n);
        }

        self_.worst_evalue = heap_array[1].best_evalue;
        if COMPO_INTENSE_DEBUG {
            assert!(s_compo_heap_is_valid(heap_array, 1, self_.n));
        }
    }

    0
}

pub fn blast_compo_heap_filled_to_cutoff(self_: &BlastCompoHeap) -> bool {
    self_.n >= self_.heap_threshold && self_.worst_evalue <= self_.ecutoff
}

pub fn blast_compo_heap_initialize(
    self_: &mut BlastCompoHeap,
    heap_threshold: i32,
    ecutoff: f64,
) -> i32 {
    self_.n = 0;
    self_.heap_threshold = heap_threshold;
    self_.ecutoff = ecutoff;
    self_.heap_array = None;
    self_.capacity = HEAP_INITIAL_CAPACITY.min(heap_threshold);
    self_.worst_evalue = 0.0;
    self_.array = Some(vec![
        BlastCompoHeapRecord {
            best_evalue: 0.0,
            best_score: 0,
            subject_index: 0,
            these_alignments: std::ptr::null_mut(),
        };
        (self_.capacity + 1).max(0) as usize
    ]);

    0
}

pub fn blast_compo_heap_release(self_: &mut BlastCompoHeap) {
    self_.heap_array = None;
    self_.array = None;
    self_.n = 0;
    self_.capacity = 0;
    self_.heap_threshold = 0;
}

pub fn blast_compo_heap_pop(self_: &mut BlastCompoHeap) -> *mut c_void {
    let mut results = std::ptr::null_mut();

    s_convert_to_heap(self_);
    if self_.n > 0 {
        let heap_array = self_
            .heap_array
            .as_mut()
            .expect("heap array is initialized");
        let first = 1usize;
        let last = self_.n as usize;

        results = heap_array[first].these_alignments;
        self_.n -= 1;
        if self_.n > 0 {
            heap_array[first] = heap_array[last];
            s_compo_heapify_down(heap_array, 1, self_.n);
        }
    }

    if COMPO_INTENSE_DEBUG {
        if let Some(heap_array) = self_.heap_array.as_ref() {
            assert!(s_compo_heap_is_valid(heap_array, 1, self_.n));
        }
    }

    results
}
