package main

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/csv"
	"fmt"
	"log"
	"math/rand"
	"os"
)

const SamplesFile = "1000genome-samples.csv"

func selectSamples(num, total int, seedPhrase string) []string {
	//	given the seed, generate N cryptographically secure pseudo-random numbers modulo M
	h := sha256.New()
	seed := int64(binary.BigEndian.Uint64(h.Sum([]byte(seedPhrase))[:8]))
	fmt.Println(seed)
	rnd := rand.New(rand.NewSource(seed))
	selection := make([]string, num)
	all := allSamples()
	for i := 0; i < num; i++ {
		selection[i] = all[rnd.Intn(total)]
	}
	return selection
}

func allSamples() []string {
	f, err := os.Open(SamplesFile)
	if err != nil {
		log.Fatalf("Error opening the samples file: %v", err)
	}
	defer f.Close()
	reader := csv.NewReader(f)
	samples, err := reader.Read()
	if err != nil {
		log.Fatalf("Error reading the samples file: %v", err)
	}
	return samples
}
