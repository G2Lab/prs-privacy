package solver

import (
	"bufio"
	"encoding/csv"
	"io"
	"log"
	"os"
	"strings"

	"github.com/nikirill/prs/data"
)

func All1000GenomesSamples() []string {
	f, err := os.Open(data.GenomesSamplesFile)
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
		"HG00124", "HG00501", "HG00635", "HG00733", "HG01983", "HG02024", "HG02046", "HG00702",
		"HG02363", "HG02372", "HG02377", "HG02381", "HG02387", "HG02388", "HG03715", "HG03948",
		"NA19240", "NA19311", "NA19313", "NA19660", "NA19675", "NA19685", "NA19985", "NA20322",
		"NA20336", "NA20341", "NA20344", "NA20526", "NA20871", "NA20893", "NA20898",
	}
}

func All1000GenomesAndRelativeSamples() []string {
	return append(All1000GenomesSamples(), AllRelativeSamples()...)
}

func AllUKBiobankSamples() []string {
	f, err := os.Open(data.UKBBSamplesFile)
	if err != nil {
		log.Fatalf("Error opening the UKBB samples file: %v", err)
	}
	defer f.Close()
	var individuals []string
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		individuals = append(individuals, scanner.Text())
	}

	if err := scanner.Err(); err != nil {
		log.Fatalf("Error reading the UKBB samples file: %v", err)
	}
	return individuals
}

func ReadRelatedIndividuals() map[string][]string {
	file, err := os.Open("info/related_individuals.txt")
	if err != nil {
		log.Fatalf("Error opening related individuals file: %v", err)
		return nil
	}
	defer file.Close()
	reader := csv.NewReader(file)
	reader.Comma = '\t'

	header, err := reader.Read()
	if err != nil {
		log.Fatalf("Error reading header: %v\n", err)
		return nil
	}
	sampleColumn, relativesColumn := -1, -1
	for i, field := range header {
		if field == "Sample" {
			sampleColumn = i
		}
		if field == "Reason for exclusion" {
			relativesColumn = i
		}
	}

	related := make(map[string][]string)
	for {
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Printf("Error reading record: %v\n", err)
			continue
		}
		split := strings.Split(record[relativesColumn], ":")
		related[record[sampleColumn]] = strings.Split(split[len(split)-1], ",")
	}
	return related
}
