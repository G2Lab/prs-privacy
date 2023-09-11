package main

import (
	"sync"
)

type MutexMap struct {
	sync.Mutex
	m map[string][]uint8
}

func NewMutexMap(regularMap map[string][]uint8) *MutexMap {
	return &MutexMap{
		m: regularMap,
	}
}

func (m *MutexMap) Get(key string) ([]uint8, bool) {
	m.Lock()
	defer m.Unlock()
	value, exists := m.m[key]
	return value, exists
}

func (m *MutexMap) Put(key string, value []uint8) {
	m.Lock()
	defer m.Unlock()
	if _, exists := m.m[key]; !exists {
		m.m[key] = make([]uint8, len(value))
		copy(m.m[key], value)
		//fmt.Printf("Added %s\n", key)
	}
}

func (m *MutexMap) RetrieveMapOnly() map[string][]uint8 {
	return m.m
}
