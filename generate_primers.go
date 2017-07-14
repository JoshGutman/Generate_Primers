package main

import "fmt"
import "os"
import "strings"
import "bufio"
import "strconv"
import "sort"


//TODO:
// > generatePrimers function
//     > Try to create the best combination of forward and reverse sequences
//     > Score the final primers (?)
// > Add the ability to allow the user to ignore a certain amount of outliers when looking for degeneracies.  For example, if there is only one
//   genome that has a degeneracy in position 1, then the user can skip it if an option is given (--ignore 1(?))
// > Command line args/options

type Sequence struct {
  name string
  seq string
  isForward bool
  val int
}


// Generates the primers
func generatePrimers(seqFile, fastaFile string) {
  degens := findDegens(seqFile, fastaFile)

  primers := make(map[string]string)

  // For each sequence
  for key, _ := range degens {
    primers[key] = ""

    // For each base in the current sequence
    for i := 0; i < len(degens[key]); i++ {

      // If the current sequence is reverse, get the code of the reverse complement
      // of the base and insert it at the beginning of the primer
      if key[len(key)-7:] == "reverse" {
        primers[key] = getCode(reverseComplement(degens[key][i])) + primers[key]

      // If the current sequence is forward, get the code of the base and put
      // it at the end of the primer
      } else {
        primers[key] += getCode(degens[key][i])
      }
    }
  }
  
  for key, val := range primers {
    fmt.Println(key + ":\t" + val)
  }
}


// Returns ambiguity code for a string of bases. The string of bases must be in
// alphabetical order. "AGT" will work but "ATG" will not.
func getCode(bases string) (string) {
  if len(bases) == 1 {
    return bases
  }
  switch bases {
    case "AG":
      return "R"
    case "CT":
      return "Y"
    case "CG":
      return "S"
    case "AT":
      return "W"
    case "GT":
      return "K"
    case "AC":
      return "M"
    case "CGT":
      return "B"
    case "AGT":
      return "D"
    case "ACT":
      return "H"
    case "ACG":
      return "V"
    case "ACGT":
      return "N"
    default:
      return "(ERROR)"
  }
}



func findDegens(seqFile, fastaFile string) (map[string]map[int]string) {

  // Key is sequence name, value is another map where key is the position of a base within the
  // sequence and value is a sorted string containing all bases that appeared in that position
  out := make(map[string]map[int]string)

  sequences := getSeqs(seqFile)
  genomes := getGenomes(fastaFile)


  // For each sequence
  for _, sequence := range sequences {

    // If sequence is a reverse, perform the reverse complement on it
    if !(sequence.isForward) {
      sequence.seq = reverseComplement(sequence.seq)
    }

    out[sequence.name] = make(map[int]string)

    // Initialize the first array slot of each sequence to the correct char
    for i, _ := range sequence.seq {
      out[sequence.name][i] += string(sequence.seq[i])
    }

    // For each genome
    for _, genome := range genomes {

      // If sequence is forward
      if sequence.isForward {

        // Iterate over each base in the sequence, looking for differences in the range specified by the sequence.val
        for i, _ := range sequence.seq {
          if sequence.seq[i] != genome[sequence.val + i - 1] && !(strings.Contains(out[sequence.name][i], string(genome[sequence.val + i - 1]))) {
            out[sequence.name][i] += string(genome[sequence.val + i - 1])
          }
        }

      // If the sequence is not forward
      } else {

        // Iterate over the sequence, performing a little bit of gymnastics to iterate forward over both the sequence.seq and the genome even when the for loop is iterating backward
        for i := len(sequence.seq) - 1; i >= 0; i-- {
          if sequence.seq[(len(sequence.seq)-1) - i] != genome[sequence.val - i - 1] && !(strings.Contains(out[sequence.name][(len(sequence.seq)-1) - i], string(genome[sequence.val - i - 1]))) {
            out[sequence.name][(len(sequence.seq)-1) - i] += string(genome[sequence.val - i - 1])
          }
        }
      }

    }

  }

  // Sort each string
  for key, val := range out {
    for i, s := range val {
      out[key][i] = sortString(s)
    }
  }

  return out
}


// Gets the nucleotide sequences from the .fasta file and
// returns them, each as single long string in a string array
func getGenomes(fastaFile string) ([]string) {

  var out []string

  // Open fasta file
  fi, err := os.Open(fastaFile)
  if err != nil {
    fmt.Println("Error - couldn't open .fasta file")
    fmt.Println(err)
    os.Exit(1)
  }
  scanner := bufio.NewScanner(fi)

  var temp string

  // Skip first line (Assuming it's a header)
  scanner.Scan()

  // For each line in the file
  for scanner.Scan() {
    line := scanner.Text()

    // If the line begins with '>', assume it's a header
    if line[0] == 62 {
      out = append(out, temp)
      temp = ""

    // If the line doesn't begin with '>', assume it's a seuence of nucleotides
    } else {
      temp += line
    }
  }

  return out
}


// Gets the sequences from the .seq file and returns them in an array of Sequence structs
// Assuming you have the following in the .seq file:
//    >99_forward
//    ACGT
// You will get a Sequence struct with the following fields:
//    name = "99_forward"
//    seq = "ACGT"
//    isForward = true
//    val = 99
func getSeqs(seqFile string) ([]Sequence) {

  var out []Sequence

  // Open the .seq file
  fi, err := os.Open(seqFile)
  if err != nil {
    fmt.Println("Error - couldn't open .seq file")
    os.Exit(1)
  }
  scanner := bufio.NewScanner(fi)

  // For each line in the file
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

    out = append(out, temp)
  }

  return out
}


// Returns the reverse complement of a sequence of nucleotides
// Only works with 'A', 'C', 'G', 'T' chars
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

/*
func contains(str []byte, token byte) (bool) {
  for _, item := range str {
    if item == token {
      return true
    }
  }
  return false
}
*/

func sortString(str string) (string) {
  if len(str) <= 1 {
    return str
  }
  s := strings.Split(str, "")
  sort.Strings(s)
  return strings.Join(s, "")
}

func main() {
  //fmt.Println(reverseComplement("GGCAAAAGCTATTTTCTCAA"))
  //fmt.Println(getSeqs("GCA_000008865.1_Escherichia_coli_Sakai_RIMD_0509952_Complete_Genome_primers.seqs"))
  //fmt.Println(findDegens("test_data\\GCA_000008865.1_Escherichia_coli_Sakai_RIMD_0509952_Complete_Genome_primers.seqs", "test_data\\all_concatenated_aligned.fasta"))
  /*
  for _, c := range []string{"A", "C", "G", "T", "AC", "AG", "AT", "CG", "CT", "GT", "ACG", "ACT", "AGT", "CGT", "ACGT"} {
    fmt.Println(getCode(c))
  }
  */
  generatePrimers("test_data\\GCA_000008865.1_Escherichia_coli_Sakai_RIMD_0509952_Complete_Genome_primers.seqs", "test_data\\all_concatenated_aligned.fasta")
}
