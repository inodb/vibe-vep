package main

import (
"encoding/gob"
"fmt"
"os"

"github.com/inodb/vibe-vep/internal/cache"
)

func main() {
f, err := os.Open(os.Args[1])
if err != nil {
panic(err)
}
defer f.Close()

var data map[string][]*cache.Transcript
if err := gob.NewDecoder(f).Decode(&data); err != nil {
panic(err)
}

gene := "CHEK2"
if len(os.Args) > 2 {
gene = os.Args[2]
}

for _, transcripts := range data {
for _, t := range transcripts {
if t.GeneName == gene {
fmt.Printf("%-35s msk=%-5v ensembl=%-5v biotype=%-30s start=%-12d end=%d\n",
t.ID, t.IsCanonicalMSK, t.IsCanonicalEnsembl, t.Biotype, t.Start, t.End)
}
}
}
}
