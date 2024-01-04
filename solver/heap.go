package solver

import "container/heap"

type genotype struct {
	mutations  []uint16
	likelihood float64
}

func newGenotype(mutations []uint16, likelihood float64) *genotype {
	mtt := make([]uint16, len(mutations))
	copy(mtt, mutations)
	return &genotype{mtt, likelihood}
}

func addToHeap(h *genheap, lkl float64, sol []uint16, heapSize int) bool {
	switch {
	case h.Len() < heapSize:
		heap.Push(h, genotype{sol, lkl})
		return true
	case lkl < (*h)[0].likelihood:
		heap.Pop(h)
		heap.Push(h, genotype{sol, lkl})
		return true
	default:
		return false
	}
}

type genheap []genotype

func (h genheap) Len() int { return len(h) }

func (h genheap) Less(i, j int) bool {
	return h[i].likelihood > h[j].likelihood
}
func (h genheap) Swap(i, j int) {
	h[i].mutations, h[j].mutations = h[j].mutations, h[i].mutations
	h[i].likelihood, h[j].likelihood = h[j].likelihood, h[i].likelihood
}

func (h *genheap) Push(x interface{}) {
	*h = append(*h, x.(genotype))
}

func (h *genheap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}
