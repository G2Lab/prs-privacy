package main

import (
	"fmt"
	"github.com/montanaflynn/stats"
	"github.com/nikirill/prs/pgs"
	"log"
	"os"
	"path/filepath"
	"sort"
)

func listFiles(folderPath string) ([]string, error) {
	dir, err := os.Open(folderPath)
	if err != nil {
		return nil, err
	}
	defer dir.Close()

	fileInfos, err := dir.Readdir(0)
	if err != nil {
		return nil, err
	}

	var fileNames []string
	for _, fileInfo := range fileInfos {
		if fileInfo.IsDir() || fileInfo.Name()[:3] != "PGS" {
			continue
		}
		fileNames = append(fileNames, fileInfo.Name())
	}

	// Sort the file names
	sort.Strings(fileNames)

	return fileNames, nil
}

func main() {
	var err error
	catalog := "catalog"
	fileNames, err := listFiles(catalog)
	if err != nil {
		log.Println("Error:", err)
		return
	}
	//fmt.Println(fileNames[:10])
	N := 1000
	numVariants := make([]float64, 0)
	for _, fileName := range fileNames[:N] {
		fmt.Println(fileName)
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join(catalog, fileName))
		if err != nil {
			log.Println("Error:", err)
			return
		}
		numVariants = append(numVariants, float64(p.VariantCount))
	}
	var mean, median float64
	median, err = stats.Median(numVariants)
	if err != nil {
		log.Println("Median error:", err)
	}
	mean, err = stats.Mean(numVariants)
	if err != nil {
		log.Println("Mean error:", err)
	}
	fmt.Printf("Median and mean number of SNPs per PGS: %d, %d\n", int(median), int(mean))
}
