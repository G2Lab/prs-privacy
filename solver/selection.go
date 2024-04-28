package solver

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
	all := All1000GenomesSamples()
	for i := 0; i < num; i++ {
		selection[i] = all[rnd.Intn(total)]
	}
	return selection
}

func All1000GenomesSamples() []string {
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

func AllRelativeSamples() []string {
	return []string{
		"HG00124", "HG00501", "HG00635", "HG00702", "HG00733", "HG01983", "HG02024", "HG02046",
		"HG02363", "HG02372", "HG02377", "HG02381", "HG02387", "HG02388", "HG03715", "HG03948",
		"NA19240", "NA19311", "NA19313", "NA19660", "NA19675", "NA19685", "NA19985", "NA20322",
		"NA20336", "NA20341", "NA20344", "NA20526", "NA20871", "NA20893", "NA20898",
	}
}

func All1000GenomesAndRelativeSamples() []string {
	return append(All1000GenomesSamples(), AllRelativeSamples()...)
}
