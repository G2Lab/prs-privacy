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
	catalog := "catalog"
	fileNames, err := listFiles(catalog)
	if err != nil {
		log.Println("Error:", err)
		return
	}
	//fmt.Println(fileNames[:10])
	//N := 100
	numVariants := make([]float64, 0)
	for _, fileName := range fileNames {
		fmt.Println(fileName)
		p := pgs.NewPGS()
		p.LoadCatalogFile(filepath.Join(catalog, fileName))
		numVariants = append(numVariants, float64(p.VariantCount))
	}
	median, err := stats.Median(numVariants)
	fmt.Printf("Median number of SNPs per PGS: %d\n", int(median))
}
