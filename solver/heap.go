package solver

import (
	"container/heap"
)

type genotype struct {
	mutations  []uint16
	likelihood float64
}

func newGenotype(mutations []uint16, likelihood float64) *genotype {
	mtt := make([]uint16, len(mutations))
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

func (h *genHeap) addToGenHeap(lkl float64, sol []uint16, heapSize int) bool {
	switch {
	case (*h).Len() < heapSize:
		ng := newGenotype(sol, lkl)
		heap.Push(h, *ng)
		return true
	case lkl < (*h)[0].likelihood:
		heap.Pop(h)
		ng := newGenotype(sol, lkl)
		heap.Push(h, *ng)
		return true
	default:
		return false
	}
}

func (h *genHeap) nthLikelihood(n int) float64 {
	return (*h)[(*h).Len()-1-n].likelihood
}
