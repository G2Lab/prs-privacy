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

type match struct {
	sums       []int64
	likelihood float64
}

func newMatch(sums []int64, likelihood float64) *match {
	return &match{sums, likelihood}
}

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

func (h *matchHeap) addToMatchHeap(lkl float64, sums []int64, heapSize int) bool {
	switch {
	case h.Len() < heapSize:
		m := newMatch(sums, lkl)
		heap.Push(h, *m)
		return true
	case lkl < (*h)[0].likelihood:
		heap.Pop(h)
		m := newMatch(sums, lkl)
		heap.Push(h, *m)
		return true
	default:
		return false
	}
}

type individual struct {
	sequence []uint8
	fitness  float64
}

func newIndividual(seq []uint8, fitness float64) *individual {
	sq := make([]uint8, len(seq))
	copy(sq, seq)
	return &individual{sq, fitness}
}

type individualHeap struct {
	individuals []*individual
	seen        map[string]struct{}
}

func newIndividualHeap() *individualHeap {
	ih := &individualHeap{make([]*individual, 0), make(map[string]struct{})}
	heap.Init(ih)
	return ih
}

func (h *individualHeap) Len() int { return len(h.individuals) }

func (h *individualHeap) Less(i, j int) bool {
	return h.individuals[i].fitness < h.individuals[j].fitness
}

func (h *individualHeap) Swap(i, j int) {
	//h.individuals[i], h.individuals[j] = h.individuals[j], h.individuals[i]
	h.individuals[i].sequence, h.individuals[j].sequence = h.individuals[j].sequence, h.individuals[i].sequence
	h.individuals[i].fitness, h.individuals[j].fitness = h.individuals[j].fitness, h.individuals[i].fitness
}

func (h *individualHeap) Push(x interface{}) {
	idv := x.(*individual)
	h.individuals = append(h.individuals, idv)
	h.seen[ArrayToString(idv.sequence)] = struct{}{}
}

func (h *individualHeap) Pop() interface{} {
	old := h.individuals
	n := len(old)
	x := old[n-1]
	old[n-1] = nil
	h.individuals = old[0 : n-1]
	delete(h.seen, ArrayToString(x.sequence))
	return x
}

func (h *individualHeap) PushIfUnseen(seq []uint8, fitness float64, heapSize int) {
	seqstr := ArrayToString(seq)
	if _, exists := h.seen[seqstr]; !exists {
		switch {
		case h.Len() < heapSize:
			idv := newIndividual(seq, fitness)
			heap.Push(h, idv)
		case h.Len() >= heapSize && fitness > h.individuals[0].fitness:
			idv := newIndividual(seq, fitness)
			heap.Pop(h)
			heap.Push(h, idv)
		default:
		}
	}
}
