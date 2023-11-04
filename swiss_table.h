#pragma once

#include <iostream>
#include <cassert>
#include <x86intrin.h>
#include <bitset>
#include <vector>
#include <utility>
#include <cstddef>
#include <stdexcept>

class BitMask {
public:
    BitMask() = default;
    explicit BitMask(const int& mask) : mask_{static_cast<unsigned long long>(mask)} {
    }
    class iterator {
    public:
        const int& operator*() const {
            return first_set_;
        }
        const int* operator->() const {
            return &first_set_;
        }
        iterator &operator++() {
            if (base_) {
                base_ &= ~(1 << first_set_);
                update();
            }
            return *this;
        }
        iterator operator++(int) {
            iterator to_return{*this};
            ++(*this);
            return to_return;
        }
        friend bool operator==(const iterator &lhs, const iterator &rhs) {
            return lhs.base_ == rhs.base_;
        }
        friend bool operator!=(const iterator &lhs, const iterator &rhs) {
            return lhs.base_ != rhs.base_;
        }
        iterator(const int parent) : base_{parent} {
            update();
        }
        void update() {
            first_set_ = __builtin_ffs(base_) - 1;
        }
        int base_;
        int first_set_;
    };
    int get() const {
        return static_cast<int>(mask_.to_ullong());
    }
    iterator begin() const {
        return {get()};
    }
    iterator end() const {
        return {0};
    }
private:
    std::bitset<16> mask_;
};

template<class KeyType, class ValueType, class Hash = std::hash<KeyType>>
class HashMap {
public:
    class SmartPair {
    public:
        using PairType = std::pair<const KeyType, ValueType>;

        SmartPair() = default;
        SmartPair(const KeyType &key, const ValueType &value) : pair_({key, value}) {}
        explicit SmartPair(std::pair<const KeyType, ValueType> pair) : pair_(pair) {}
        PairType &get_pair() {
            return pair_;
        }
        const PairType &get_pair() const {
            return pair_;
        }
        KeyType get_first() const {
            return pair_.first;
        }
        ValueType get_second() const {
            return pair_.second;
        }
        ValueType &get_ref_second() {
            return pair_.second;
        }
        const ValueType &get_ref_second() const {
            return pair_.second;
        }
        std::pair<KeyType &, ValueType &> get_ref() {
            return std::pair<KeyType &, ValueType &>(const_cast<KeyType &>(pair_.first), pair_.second);
        }
        SmartPair &operator=(const SmartPair &other) {
            get_ref() = other.pair_;
            return *this;
        }
    private:
        PairType pair_;
    };
    struct slot {
        SmartPair element;
        bool empty = true;
        bool deleted = false;

        slot() = default;
        slot(KeyType key, ValueType value) : element(key, value), empty(false) {}
        explicit slot(std::pair<KeyType, ValueType> other_element) : element(other_element), empty(false) {}
        slot(const slot &other) {
            element = other.element;
            empty = other.empty;
            deleted = other.deleted;
        }
        KeyType get_key() const {
            return element.get_first();
        }
        ValueType get_value() const {
            return element.get_second();
        }
        slot &operator=(const slot &other) {
            element = other.element;
            empty = other.empty;
            deleted = other.deleted;
            return *this;
        }
        bool operator==(const slot &other) {
            return element.get_pair() == other.element.get_pair() and empty == other.empty and deleted == other.deleted;
        }
    };
    class iterator {
    public:
        using IteratorType = typename std::vector<slot>::iterator;

        iterator() = default;
        iterator(IteratorType pointer, IteratorType end) :
                pointer_(pointer), end_(end) {
        }
        std::pair<const KeyType, ValueType> &operator*() const {
            return pointer_->element.get_pair();
        }
        std::pair<const KeyType, ValueType> *operator->() const {
            return &(pointer_->element.get_pair());
        }
        iterator &operator++() {
            ++pointer_;
            while (pointer_ != end_ and (pointer_->empty or pointer_->deleted)) {
                ++pointer_;
            }
            return *this;
        }
        iterator operator++(int) {
            iterator copy_iterator(*this);
            ++pointer_;
            while (pointer_ != end_ and (pointer_->empty or pointer_->deleted)) {
                ++pointer_;
            }
            return copy_iterator;
        }
        iterator operator+(size_t shift) const {
            iterator copy_iterator(*this);
            copy_iterator.pointer_ += shift;
            return copy_iterator;
        }
        size_t operator-(const iterator &other) const {
            return pointer_ - other.pointer_;
        }
        bool operator==(const iterator &other) const {
            return pointer_ == other.pointer_;
        }
        bool operator!=(const iterator &other) const {
            return pointer_ != other.pointer_;
        }
    private:
        IteratorType pointer_;
        IteratorType end_;
    };
    class const_iterator {
    public:
        using ConstIteratorType = typename std::vector<slot>::const_iterator;

        const_iterator() = default;
        const_iterator(ConstIteratorType pointer, ConstIteratorType end) : pointer_(pointer), end_(end) {
        }
        const std::pair<const KeyType, ValueType> &operator*() const {
            return pointer_->element.get_pair();
        }
        const std::pair<const KeyType, ValueType> *operator->() const {
            return &(pointer_->element.get_pair());
        }
        const_iterator &operator++() {
            ++pointer_;
            while (pointer_ != end_ and (pointer_->empty or pointer_->deleted)) {
                ++pointer_;
            }
            return *this;
        }
        const_iterator operator++(int) {
            const_iterator copy_iterator(*this);
            ++pointer_;
            while (pointer_ != end_ and (pointer_->empty or pointer_->deleted)) {
                ++pointer_;
            }
            return copy_iterator;
        }
        size_t operator-(const const_iterator &other) const {
            return pointer_ - other.pointer_;
        }
        bool operator==(const const_iterator &other) const {
            return pointer_ == other.pointer_;
        }
        bool operator!=(const const_iterator &other) const {
            return pointer_ != other.pointer_;
        }
    private:
        ConstIteratorType pointer_;
        ConstIteratorType end_;
    };

    explicit HashMap(Hash hash = Hash())
            : hash_(hash),
              hash_map_(std::vector<slot>(sizes_.front() * 16)),
              control_bytes_(std::vector<uint8_t>(sizes_.front() * 16, 0b11111111)) {}
    template<class Iterator>
    HashMap(Iterator begin, Iterator end, Hash hash = Hash())
            : hash_(hash),
              hash_map_(std::vector<slot>(sizes_.front() * 16)),
              control_bytes_(std::vector<uint8_t>(sizes_.front() * 16, 0b11111111)) {
        while (begin != end) {
            insert(*begin);
            ++begin;
        }
    }
    HashMap(std::initializer_list<std::pair<KeyType, ValueType>> seq, Hash hash = Hash())
            : hash_(hash),
              hash_map_(std::vector<slot>(sizes_.front() * 16)),
              control_bytes_(std::vector<uint8_t>(sizes_.front() * 16, 0b11111111)) {
        for (auto element : seq) {
            insert(element);
        }
    }

    size_t size() const {
        return real_size_;
    }
    bool empty() const {
        return real_size_ == 0;
    }
    Hash hash_function() const {
        return hash_;
    }

    void insert(std::pair<KeyType, ValueType> element) {
        KeyType key = element.first;

        if (find(key) == end()) {
            ++real_size_;
            ++captured_cnt_;
            size_t group = get_hash1(key) % sizes_[cur_size_];
            int ITERS = 0;
            while (true) {
                ITERS++; if (ITERS >= 1e6) assert(false);
                BitMask i = match(0b11111111, group * 16);

                if (i.begin() != i.end()) {
                    hash_map_[group * 16 + *i.begin()] = slot(element);
                    control_bytes_[group * 16 + *i.begin()] = get_hash2(key);
                    break;
                }
                i = match(0b10000000, group * 16);
                if (i.begin() != i.end()) {
                    hash_map_[group * 16 + *i.begin()] = slot(element);
                    control_bytes_[group * 16 + *i.begin()] = get_hash2(key);
                    break;
                }
                group = (group + 1) % sizes_[cur_size_];
            }
            if (static_cast<double>(captured_cnt_) / hash_map_.size() >= MAX_LOAD_FACTOR) {
                rebuild();
            }
        }
    }
    void erase(KeyType key) {
        iterator it = find(key);
        if (it != end()) {
            --real_size_;
            auto cur_begin = iterator(hash_map_.begin(), hash_map_.end());
            int group = (it - cur_begin) / 16;
            if (check_empty(group * 16)) {
                hash_map_[it - cur_begin].empty = true;
                control_bytes_[it - cur_begin] = 0b11111111;
            } else {
                hash_map_[it - cur_begin].deleted = true;
                control_bytes_[it - cur_begin] = 0b10000000;
            }
        }
    }
    iterator find(KeyType key) {
        size_t group = get_hash1(key) % sizes_[cur_size_];
        int ITERS = 0;
        while (true) {
            ITERS++; if (ITERS >= 1e6) assert(false);
            for (auto i : match(get_hash2(key), group * 16)) {
                if (key == hash_map_[group * 16 + i].get_key()) {
                    return iterator(hash_map_.begin() + group * 16 + i, hash_map_.end());
                }
            }
            if (check_empty(group * 16)) return end();
            group = (group + 1) % sizes_[cur_size_];
        }
    }
    const_iterator find(KeyType key) const {
        size_t group = get_hash1(key) % sizes_[cur_size_];
        int ITERS = 0;
        while (true) {
            ITERS++; if (ITERS >= 1e6) assert(false);
            for (auto i : match(get_hash2(key), group * 16)) {
                if (key == hash_map_[group * 16 + i].get_key()) {
                    return const_iterator(hash_map_.cbegin() + group * 16 + i, hash_map_.cend());
                }
            }
            if (check_empty(group * 16)) return end();
            group = (group + 1) % sizes_[cur_size_];
        }
    }
    ValueType &operator[](KeyType key) {
        insert({key, ValueType()});
        iterator it = find(key);
        return it->second;
    }
    const ValueType &at(KeyType key) const {
        const_iterator it = find(key);
        if (it == end()) {
            throw std::out_of_range("There is no such key in hash map");
        }
        return it->second;
    }
    void clear() {
        cur_size_ = 0;
        real_size_ = 0;
        captured_cnt_ = 0;
        hash_map_.clear();
        hash_map_.resize(sizes_.front() * 16);
        control_bytes_.clear();
        control_bytes_.resize(sizes_.front() * 16, 0b11111111);
    }

    iterator begin() {
        for (size_t i = 0; i < hash_map_.size(); ++i) {
            if (control_bytes_[i] != 0b11111111 and control_bytes_[i] != 0b10000000) {
                return iterator(hash_map_.begin() + i, hash_map_.end());
            }
        }
        return iterator(hash_map_.end(), hash_map_.end());
    }
    iterator end() {
        return iterator(hash_map_.end(), hash_map_.end());
    }
    const_iterator begin() const {
        for (size_t i = 0; i < hash_map_.size(); ++i) {
            if (control_bytes_[i] != 0b11111111 and control_bytes_[i] != 0b10000000) {
                return const_iterator(hash_map_.cbegin() + i, hash_map_.cend());
            }
        }
        return const_iterator(hash_map_.cend(), hash_map_.cend());
    }
    const_iterator end() const {
        return const_iterator(hash_map_.cend(), hash_map_.cend());
    }
private:
    constexpr const static double MAX_LOAD_FACTOR = 0.6;
    std::vector<int> sizes_ =
            {5, 11, 23, 47, 97, 197, 397, 797, 1597, 3203, 6421, 12853, 25717, 51437, 102877, 205759, 411527, 823117,
             1646237, 3292489, 6584983};
    size_t cur_size_ = 0;
    size_t real_size_ = 0;
    size_t captured_cnt_ = 0;
    Hash hash_;
    std::vector<slot> hash_map_;
    std::vector<uint8_t> control_bytes_;

    BitMask match(uint8_t hash, size_t index) const {
        auto distributed_hash = _mm_set1_epi8(hash);
        auto metadata = _mm_loadu_si128((__m128i *) &control_bytes_[index]);
        return BitMask(_mm_movemask_epi8(_mm_cmpeq_epi8(distributed_hash, metadata)));
    }
    bool check_empty(size_t index) const {
        uint8_t hash = 0b11111111;
        auto distributed_hash = _mm_set1_epi8(hash);
        auto metadata = _mm_loadu_si128((__m128i *) &control_bytes_[index]);
        return _mm_movemask_epi8(_mm_cmpeq_epi8(distributed_hash, metadata)) != 0;
    }

    size_t get_hash1(const KeyType &key) const {
        return hash_(key) >> 7;
    }

    size_t get_hash2(const KeyType &key) const {
        return hash_(key) & 127;
    }

    void rebuild() {
        std::vector<slot> copy_hash_map = hash_map_;
        std::vector<uint8_t> copy_control_bytes = control_bytes_;

        real_size_ = 0;
        captured_cnt_ = 0;
        ++cur_size_;
        hash_map_.clear();
        hash_map_.resize(sizes_[cur_size_] * 16);
        control_bytes_.clear();
        control_bytes_.resize(sizes_[cur_size_] * 16, 0b11111111);
        for (size_t i = 0; i < copy_hash_map.size(); ++i) {
            if (copy_control_bytes[i] != 0b11111111 and copy_control_bytes[i] != 0b10000000) {
                insert({copy_hash_map[i].get_key(), copy_hash_map[i].get_value()});
            }
        }
    }
};
