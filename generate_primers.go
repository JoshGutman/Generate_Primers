package main

import "fmt"
import "os"
import "strings"
import "bufio"
import "strconv"


//TODO:
// > Sort each array (the value of the nested map) in findDegens
// > Return the nested maps in findDegens
// > generate_primers function -- take in the nested maps and generate a few of the best primers
//     > Best = try to pick the reverse and forward reads with the fewest degeneracies
//     > Once the best combinations are picked, convert them to the format that includes the arbitrary bases (e.g. A OR T OR C = H or whatever)
//     > Don't forget to reverse-Complement the reverse sequences again before doing the above step
//     > Score the final primers (?)
// > Add the ability to allow the user to ignore a certain amount of outliers when looking for degeneracies.  For example, if there is only one
//   genome that has a degeneracy in position 1, then the user can skip it if an option is given (--fuzz 1(?))
// > Command line args/options 


/*

Takes in target_primers.seqs and all_concatenated_alligned.fasta

For each sequence in target_primers.seq:
  > Calculate degeneracies
    >> Be sure to keep track of which bases are in the degeneracies
    >> Be sure to do reverse complement of reverse sequences
    >> Store in map
  > Calculate primer

Choose the best forward and reverse sequence (make sure 100 < |XX -YY| + 1 < 250)


*/

type Sequence struct {
  name string
  seq string
  isForward bool
  val int
  length int
}


func findDegens(seqFile, fastaFile string) {

  // Key is sequence name, value is another map where key is the position of a base within the sequence and value is
  //  a sorted byte containing the numeric representation (type=byte) of all bases that appeared in that position
  out := make(map[string]map[int][]byte)

  sequences := getSeqs(seqFile)
  genomes := getGenomes(fastaFile)


  // For each sequence
  for _, sequence := range sequences {

    // If sequence is a reverse, perform the reverse complement on it
    if !(sequence.isForward) {
      sequence.seq = reverseComplement(sequence.seq)
    }

    out[sequence.name] = make(map[int][]byte)

    // Initialize the first array slot of each sequence to the correct char
    for i, _ := range sequence.seq {
      out[sequence.name][i] = append(out[sequence.name][i], sequence.seq[i])
    }

    // For each genome
    for _, genome := range genomes {

      // If sequence is forward
      if sequence.isForward {

        // Iterate over each base in the sequence, looking for differences in the range specified by the sequence.val
        for i, _ := range sequence.seq {
          if sequence.seq[i] != genome[sequence.val + i - 1] && !(contains(out[sequence.name][i], genome[sequence.val + i - 1])) {
            out[sequence.name][i] = append(out[sequence.name][i], genome[sequence.val + i - 1])
          }
        }

      // If the sequence is not forward
      } else {

        // Iterate over the sequence, performing a little bit of gymnastics to iterate forward over both the sequence.seq and the genome even when the for loop is iterating backward
        for i := len(sequence.seq) - 1; i >= 0; i-- {
          if sequence.seq[(len(sequence.seq)-1) - i] != genome[sequence.val - i - 1] && !(contains(out[sequence.name][(len(sequence.seq)-1) - i], genome[sequence.val - i - 1])) {
            out[sequence.name][(len(sequence.seq)-1) - i] = append(out[sequence.name][(len(sequence.seq)-1) - i], genome[sequence.val - i - 1])
          }
        }
      }

    }

  }

  for key, val := range out {
    fmt.Println(key)
    for _, i := range val {
      fmt.Println(i)
    }
    fmt.Println()
  }
}


func getGenomes(fastaFile string) ([]string) {

  var out []string

  fi, err := os.Open(fastaFile)
  if err != nil {
    fmt.Println("Error - couldn't open .fasta file")
    fmt.Println(err)
    os.Exit(1)
  }
  scanner := bufio.NewScanner(fi)

  var temp string

  scanner.Scan()

  for scanner.Scan() {
    line := scanner.Text()
    if line[0] == 62 {
      out = append(out, temp)
      temp = ""
    } else {
      temp += line
    }
  }

  return out
}


func getSeqs(seqFile string) ([]Sequence) {

  var out []Sequence

  fi, err := os.Open(seqFile)
  if err != nil {
    fmt.Println("Error - couldn't open .seq file")
    os.Exit(1)
  }
  scanner := bufio.NewScanner(fi)

  for scanner.Scan() {

    var temp Sequence

    // Get name
    line := scanner.Text()[1:]
    temp.name = line

    // Get value
    split_line := strings.Split(line, "_")
    temp.val, _ = strconv.Atoi(split_line[0])

    // Get isForward
    if split_line[1] == "forward" {
      temp.isForward = true
    } else {
      temp.isForward = false
    }

    // Get sequence
    scanner.Scan()
    temp.seq = scanner.Text()

    // Get length of sequence
    temp.length = len(temp.seq)

    out = append(out, temp)
  }

  return out
}



func reverseComplement(sequence string) (out string) {
  for i := len(sequence)-1; i >= 0; i-- {

    switch sequence[i] {

    case 65:
      out += "T"
      break
    case 84:
      out += "A"
      break
    case 71:
      out += "C"
      break
    case 67:
      out += "G"
      break
    default:
      fmt.Println("Error -- Encountered non-ATGC char in sequence")
    }

  }
  return
}


func contains(str []byte, token byte) (bool) {
  for _, item := range str {
    if item == token {
      return true
    }
  }
  return false
}

func main() {
  //fmt.Println(reverseComplement("GGCAAAAGCTATTTTCTCAA"))
  //fmt.Println(getSeqs("GCA_000008865.1_Escherichia_coli_Sakai_RIMD_0509952_Complete_Genome_primers.seqs"))
  //fmt.Println(contains([]int {1, 2, 3}, 4))
  findDegens("GCA_000008865.1_Escherichia_coli_Sakai_RIMD_0509952_Complete_Genome_primers.seqs", "all_concatenated_aligned.fasta")
}
