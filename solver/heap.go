package solver

import "container/heap"

type genotype struct {
	mutations  []uint8
	likelihood float32
}

func newGenotype(mutations []uint8, likelihood float32) *genotype {
	mtt := make([]uint8, len(mutations))
	copy(mtt, mutations)
	return &genotype{mtt, likelihood}
}

type genHeap []genotype

func newGenHeap() *genHeap {
	h := make(genHeap, 0)
	heap.Init(&h)
	return &h
}

func (h *genHeap) Len() int { return len(*h) }

func (h *genHeap) Less(i, j int) bool {
	return (*h)[i].likelihood > (*h)[j].likelihood
}
func (h *genHeap) Swap(i, j int) {
	(*h)[i].mutations, (*h)[j].mutations = (*h)[j].mutations, (*h)[i].mutations
	(*h)[i].likelihood, (*h)[j].likelihood = (*h)[j].likelihood, (*h)[i].likelihood
}

func (h *genHeap) Push(x interface{}) {
	*h = append(*h, x.(genotype))
}

func (h *genHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

func (h *genHeap) addToGenHeap(lkl float32, sol []uint8, heapSize int) bool {
	if (*h).Len() < heapSize {
		ng := newGenotype(sol, lkl)
		heap.Push(h, *ng)
		return true
	}
	if lkl < (*h)[0].likelihood {
		heap.Pop(h)
		ng := newGenotype(sol, lkl)
		heap.Push(h, *ng)
		return true
	}
	return false
}

type match struct {
	sums       []int64
	likelihood float32
}

func newMatch(sums []int64, likelihood float32) *match {
	return &match{sums, likelihood}
}

// Implement as maxheap, i.e., the first element is the largest
type matchHeap []match

func newMatchHeap() *matchHeap {
	h := make(matchHeap, 0)
	heap.Init(&h)
	return &h
}

func (h *matchHeap) Len() int { return len(*h) }

func (h *matchHeap) Less(i, j int) bool {
	return (*h)[i].likelihood > (*h)[j].likelihood
}
func (h *matchHeap) Swap(i, j int) {
	(*h)[i].sums, (*h)[j].sums = (*h)[j].sums, (*h)[i].sums
	(*h)[i].likelihood, (*h)[j].likelihood = (*h)[j].likelihood, (*h)[i].likelihood
}

func (h *matchHeap) Push(x interface{}) {
	*h = append(*h, x.(match))
}

func (h *matchHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

func (h *matchHeap) addToMatchHeap(lkl float32, sums []int64, heapSize int) bool {
	if h.Len() < heapSize {
		m := newMatch(sums, lkl)
		heap.Push(h, *m)
		return true
	}
	if lkl < (*h)[0].likelihood {
		heap.Pop(h)
		m := newMatch(sums, lkl)
		heap.Push(h, *m)
		return true
	}
	return false
}
