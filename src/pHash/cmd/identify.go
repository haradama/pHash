package cmd

import (
	"bytes"
	"fmt"
	"html/template"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"runtime"
	"strings"
	"sync"

	"github.com/OneOfOne/xxhash"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
)

func init() {
	RootCmd.AddCommand(identifyCmd)
	identifyCmd.Flags().StringVarP(&o.optIn, "in", "i", "", "Input FASTA file")
	identifyCmd.Flags().StringVarP(&o.optDB, "db", "d", "", "Database")
	identifyCmd.Flags().IntVarP(&o.optThreshold, "threshold", "t", 10, "Threshold of probability")
}

var identifyCmd = &cobra.Command{
	Use:   "identify",
	Short: "Identifier of plasmid",
	Long:  "Identifier of plasmid using database",
	Run: func(cmd *cobra.Command, args []string) {

		if o.optIn == "" || o.optDB == "" {
			fmt.Println("--in and --db are required")
			cmd.Help()
			os.Exit(0)
		}

		inFile := o.optIn
		db := o.optDB
		threshold := float32(o.optThreshold) * 0.01
		outFile := "pHash.log.txt"

		cpus := runtime.NumCPU()
		runtime.GOMAXPROCS(cpus)

		var in *fasta.Reader
		if inFile == "" {
			in = fasta.NewReader(os.Stdin, linear.NewSeq("", nil, alphabet.DNA))
		} else if f, err := os.Open(inFile); err != nil {
			fmt.Fprintf(os.Stderr, "Error: %v.", err)
			os.Exit(1)
		} else {
			in = fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA))
			defer f.Close()
		}

		ambiguousDnaComplement := strings.NewReplacer(
			"A", "T",
			"C", "G",
			"G", "C",
			"T", "A",
			"M", "K",
			"R", "Y",
			"Y", "R",
			"K", "M",
			"V", "B",
			"H", "D",
			"D", "H",
			"B", "V")

		type Row struct {
			AccID   string
			Seq     string
			Length  int
			Link    template.HTML
			Phylum  string
			Jaccard float32
		}

		type table []Row
		var tb table

		binary, err := ioutil.ReadFile(db)
		if err != nil {
			if err != io.EOF {
				fmt.Println(err)
				os.Exit(1)
			}
			return
		}

		k := messagePackDecoding(&binary).Kmer
		sketchSize := messagePackDecoding(&binary).SketchSize
		plasmidsRecords := messagePackDecoding(&binary).Plasmid

		var wg sync.WaitGroup
		mutex := new(sync.Mutex)

		fwfasta, err := os.Create("pHash_plasmids.fna")
		if err != nil {
			if err != io.EOF {
				fmt.Println(err)
				os.Exit(1)
			}
			return
		}
		fastaw := fasta.NewWriter(fwfasta, 60)
		plasmidSeq := linear.NewSeq("", nil, alphabet.DNA)

		fw, err := os.Create(outFile)
		if err != nil {
			if err != io.EOF {
				fmt.Println(err)
				os.Exit(1)
			}
			return
		}
		defer fw.Close()

		line := "AccId\tSimilarPlasmidAccId\tSimilarity\n"
		fw.Write(([]byte)(line))

		for i := 0; i < 1000; i++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				for {
					mutex.Lock()
					s, err := in.Read()
					mutex.Unlock()

					if err != nil {
						if err != io.EOF {
							fmt.Println(err)
							os.Exit(1)
						}
						break
					}

					read := s.Slice()
					kmerNum := read.Len() - (k - 1)
					kmerMap := make(map[string]struct{}, kmerNum)
					var seq string

					for i := 0; i < kmerNum; i++ {
						seq = read.Slice(i, i+k).(alphabet.Letters).String()
						if strings.Index(seq, "N") == -1 {
							reverseComplement := ambiguousDnaComplement.Replace(rev(&seq))

							h1 := xxhash.ChecksumString32S(seq, 0)
							h2 := xxhash.ChecksumString32S(reverseComplement, 0)

							var canonicalKmer string
							if h1 > h2 {
								canonicalKmer = seq
							} else {
								canonicalKmer = reverseComplement
							}
							kmerMap[canonicalKmer] = struct{}{}
						}
					}

					kmerList := make([][]byte, len(kmerMap))
					for key := range kmerMap {
						kmerList = append(kmerList, []byte(key))
					}

					minHashValues := make([]*uint64, sketchSize)

					for i := uint64(0); i < sketchSize; i++ {
						minValue := uint64((1 << 64) - 1)
						for _, kmer := range kmerList {
							xxhv := xxhash.Checksum64S(kmer, i)
							if minValue > xxhv {
								minValue = xxhv
							}
						}
						minHashValues[i] = &minValue
					}

					var (
						bestHitKeys  []string
						bestHitValue float32
						similarity   float32
					)

					for _, record := range plasmidsRecords {
						refKey := record.AccID
						refValue := record.PlasmidMinHashValue
						similarity = calcSimilarity(minHashValues, &refValue, sketchSize)

						if similarity > bestHitValue {
							bestHitKeys = []string{refKey}
							bestHitValue = similarity
						} else if similarity == bestHitValue {
							bestHitKeys = append(bestHitKeys, refKey)
						}
					}

					mutex.Lock()
					seqSymbol := fmt.Sprintf("%s...%s", seq[:3], seq[len(seq)-3:])

					var plasmidURL string
					if len(bestHitKeys) == 1 {
						plasmidURL = fmt.Sprintf("<a href=\"https://www.ncbi.nlm.nih.gov/nuccore/%s\" target=\"_blank\">%s</a>", bestHitKeys[0], bestHitKeys[0])

						line = fmt.Sprintf("%s\t%s\t%f\n", s.Name(), bestHitKeys[0], bestHitValue)
						fw.Write(([]byte)(line))

					} else if len(bestHitKeys) > 1 {
						for _, bestHitKey := range bestHitKeys {
							plasmidURLs := []string{}
							plasmidURLs = append(plasmidURLs, fmt.Sprintf("<a href=\"https://www.ncbi.nlm.nih.gov/nuccore/%s\" target=\"_blank\">%s</a>", bestHitKey, bestHitKeys[0]))
							plasmidURL = strings.Join(plasmidURLs, "<br>")

							line = fmt.Sprintf("%s\t%s\t%f\n", s.Name(), bestHitKey, bestHitValue)
							fw.Write(([]byte)(line))
						}
					}

					if bestHitValue >= threshold {
						tb = append(tb, Row{AccID: s.Name(), Seq: seqSymbol, Length: s.Len(), Link: template.HTML(plasmidURL), Phylum: "Proteobacteria", Jaccard: bestHitValue})

						plasmidSeq.ID = s.Name()
						plasmidSeq.Seq = s.Slice().(alphabet.Letters)
						plasmidSeq.Desc = fmt.Sprintf("Similar to %s (%f)", strings.Join(bestHitKeys, ":"), bestHitValue)

						_, err := fastaw.Write(plasmidSeq)
						if err != nil {
							log.Printf("Failed to write: %s", err)
						}
					}

					mutex.Unlock()
				}
			}()
		}
		wg.Wait()

		if err := os.MkdirAll("report/assets", 0777); err != nil {
			if err != io.EOF {
				fmt.Println(err)
				os.Exit(1)
			}
		}
		copyFile("/assets/assets/pHash_logo.svg", "report/assets/pHash_logo.svg")
		copyFile("/assets/assets/bootstrap.bundle.min.js", "report/assets/bootstrap.bundle.min.js")
		copyFile("/assets/assets/bootstrap.min.css", "report/assets/bootstrap.min.css")
		copyFile("/assets/assets/bootstrap.min.js", "report/assets/bootstrap.min.js")

		f, err := Assets.Open("/assets/template.html.tpl")
		if err != nil {
			if err != io.EOF {
				fmt.Println(err)
				os.Exit(1)
			}
		}
		defer f.Close()

		contents, _ := ioutil.ReadAll(f)

		t, err := template.New("").Parse(string(contents))
		if err != nil {
			if err != io.EOF {
				fmt.Println(err)
				os.Exit(1)
			}
		}

		report, err := os.Create("report/index.html")
		if err != nil {
			if err != io.EOF {
				fmt.Println(err)
				os.Exit(1)
			}
		}
		defer report.Close()

		buff := new(bytes.Buffer)
		fwt := io.Writer(buff)

		if err := t.Execute(fwt, tb); err != nil {
			log.Fatal(err)
		}
		report.Write(buff.Bytes())

	},
}

func calcSimilarity(array1 []*uint64, array2 *[]uint64, size uint64) float32 {
	var match float64
	for i := 0; i < int(size); i++ {
		if *array1[i] == (*array2)[i] {
			match++
		}
	}
	var bitNum float64 = 64
	similarity := (match / float64(size)) - math.Pow(2, -bitNum)/(1-math.Pow(2, -bitNum))

	return float32(similarity)
}

func copyFile(srcName string, dstName string) {
	src, err := Assets.Open(srcName)
	if err != nil {
		if err != io.EOF {
			fmt.Println(err)
			os.Exit(1)
		}
	}
	defer src.Close()

	dst, err := os.Create(dstName)
	if err != nil {
		if err != io.EOF {
			fmt.Println(err)
			os.Exit(1)
		}
	}
	defer dst.Close()

	_, err = io.Copy(dst, src)
	if err != nil {
		if err != io.EOF {
			fmt.Println(err)
			os.Exit(1)
		}
	}
}
