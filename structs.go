package main

import (
	"sync"
)

type MutexMap struct {
	sync.Mutex
	m map[string][]int
}

func NewMutexMap(regularMap map[string][]int) *MutexMap {
	return &MutexMap{
		m: regularMap,
	}
}

func (m *MutexMap) Get(key string) ([]int, bool) {
	m.Lock()
	defer m.Unlock()
	value, exists := m.m[key]
	return value, exists
}

func (m *MutexMap) Put(key string, value []int) {
	m.Lock()
	defer m.Unlock()
	if _, exists := m.m[key]; !exists {
		m.m[key] = make([]int, len(value))
		copy(m.m[key], value)
		//fmt.Printf("Added %s\n", key)
	}
}

func (m *MutexMap) RetrieveMapOnly() map[string][]int {
	return m.m
}
