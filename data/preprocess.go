package main

import (
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"os"
	"strings"

	"github.com/nikirill/prs/solver"
	"github.com/nikirill/prs/tools"
)

func readSamplePopulations() {
	fileName := "data/igsr_samples.tsv"

	// Open the TSV inputFile
	inputFile, err := os.Open(fileName)
	if err != nil {
		fmt.Println("Error:", err)
		return
	}
	defer inputFile.Close()

	// Create a CSV reader with tab as the delimiter
	reader := csv.NewReader(inputFile)
	reader.Comma = '\t'

	// Read the header to get column names
	header, err := reader.Read()
	if err != nil {
		fmt.Println("Error reading header:", err)
		return
	}

	sampleIdx := indexOf(header, "Sample name")
	superPopIdx := indexOf(header, "Superpopulation code")

	if sampleIdx == -1 || superPopIdx == -1 {
		fmt.Println("Required fields not found in the header.")
		return
	}

	fullData := make(map[string]string)
	var record []string
	// Read the remaining records and extract fields
	for {
		record, err = reader.Read()
		if err != nil {
			break // End of inputFile
		}
		fullData[record[sampleIdx]] = record[superPopIdx]
	}

	//// Check that we have all the sample info for the 1000 Genomes dataset
	//sampleFile, err := os.Open("data/1000genome-samples.csv")
	//if err != nil {
	//	fmt.Println("Error:", err)
	//	return
	//}
	//defer sampleFile.Close()
	//
	//// Create a CSV reader with tab as the delimiter
	//reader = csv.NewReader(sampleFile)
	//
	//// Read the header to get column names
	//sampleList, err := reader.Read()
	//if err != nil {
	//	fmt.Println("Error reading sample list:", err)
	//	return
	//}
	//reducedData := make(map[string]string)
	//for _, sample := range sampleList {
	//	if _, exists := fullData[sample]; !exists {
	//		fmt.Printf("Sample %s not found in the superpopulation fullData\n", sample)
	//		continue
	//	}
	//	reducedData[sample] = fullData[sample]
	//}

	// Save the results to a JSON file
	outputFile, err := os.Create("data/superpopulations.json")
	if err != nil {
		fmt.Println("Error creating JSON outputFile:", err)
		return
	}
	defer outputFile.Close()

	// Encode the fullData as JSON and write to the outputFile
	encoder := json.NewEncoder(outputFile)
	err = encoder.Encode(fullData)
	if err != nil {
		fmt.Println("Error encoding fullData to JSON:", err)
	}
}

// indexOf finds the index of a value in a slice, or returns -1 if not found
func indexOf(slice []string, value string) int {
	for i, v := range slice {
		if v == value {
			return i
		}
	}
	return -1
}

func saveSamplesWithoutRelatives() {
	relations := solver.ReadRelatedIndividuals()
	haveExternalRelatives := make(map[string]struct{})
	for _, relatives := range relations {
		for _, relative := range relatives {
			haveExternalRelatives[relative] = struct{}{}
		}
	}
	allSamples := solver.All1000GenomesSamples()
	ancestry := tools.LoadAncestry()
	withoutRelatives := make(map[string][]string)
	var exists bool
	var ppl string
	for _, sample := range allSamples {
		if _, exists = haveExternalRelatives[sample]; exists {
			fmt.Printf("Sample %s is excluded due to external relatives\n", sample)
			continue
		}
		ppl = ancestry[sample]
		if strings.Contains(ppl, ",") {
			fmt.Printf("Sample %s is excluded due to multi-ancestry %s\n", sample, ppl)
			continue
		}
		if _, exists = withoutRelatives[ppl]; !exists {
			withoutRelatives[ppl] = make([]string, 0)
		}
		withoutRelatives[ppl] = append(withoutRelatives[ppl], sample)
	}
	for population := range withoutRelatives {
		file, err := os.OpenFile(fmt.Sprintf("data/1000g_%s_no_relatives.txt", population), os.O_CREATE|os.O_WRONLY, 0644)
		if err != nil {
			log.Fatalf("Cannot create file: %v", err)
		}
		defer file.Close()
		for i, sample := range withoutRelatives[population] {
			if i != 0 {
				_, err = file.WriteString("\n")
				if err != nil {
					log.Fatalf("Cannot write to file: %v", err)
				}
			}
			_, err = file.WriteString(sample)
			if err != nil {
				log.Fatalf("Cannot write to file: %v", err)
			}
		}
	}
}

func main() {
	readSP := flag.Bool("populations", false, "read sample populations")
	relatives := flag.Bool("relatives", false, "save samples without relatives")
	flag.Parse()

	if *readSP {
		readSamplePopulations()
	}
	if *relatives {
		saveSamplesWithoutRelatives()
	}
}
