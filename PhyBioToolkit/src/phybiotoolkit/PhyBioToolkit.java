/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package phybiotoolkit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import joey.SeqPaint.SeqPaint;
import joey.Structure.ParseStructure;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.ProteinSequence;

/**
 *
 * @author zhongxf
 */
public class PhyBioToolkit {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IllegalSymbolException, FileNotFoundException, IOException 
    {
//        testSeqPaint();
        //JAL-1266 extract RMSDs from Jmol superposition of alignment
        testParsStructure(); 
    }

    private static void testSeqPaint() throws IllegalSymbolException 
    {
        Sequence[] nrSeqs = new Sequence[3];
        ProteinSequence[] aaSeqs = new ProteinSequence[3];
        try {
////            String tmpString1 = "GAAATCGATCGATAGCTTTTTTTTTTTACGATA-GACTAGCATTCCGACGATA-GACTAGCATTCCC";
////            String tmpString2 = "AAAATCGATC-ATAGC----TT----TACGATACGACTAGCATTCCGAC--TA-GACTAGCATTCC-";
            //String tmpString3 = "GAAAT--ATC-ATAGC----------TACGATACGACTAGCATTCCGAC--TA--ACTAGG----CC";
            String tmpString1 = "GAGATAGGGTGGCAGATGTAATTGAAAGTTCCATAGGAGATAGCGTG-GCAGAGCCCTCACTCACGCTCTACCAGCACCCACAGGCCAG-ACACACAGGTGAGCAGTCATCGACTGGA-ACAGGCAAGGTTCCAGCACTCCAAGCTGCTGAAATTGGAGCATCATCAAATGCTAGTGACGAGAG-ATGATTGAGACACGCTGTGT-CTTAACTCGCACAGCACAGCTGAGAC-ACTCTTGATAGTTTCTTCAGCAGAGCGGGATTAGTTGGAGAG-TAGATCTCCCTCTTGAAGGCACAACTAACCCAAATGGTTATGCCAACTGGGACATAGATATAACAGGTTACGCGCAAATGCGTAGAAAGGTAGAGCT-TTCACCTACATGCGCTTTGATGCAGAGTTCACTTTTGTTGCGTGCACACCCACCGGG--AGTTGTCCCACAATTGCTCCAATATATGTTTGTGCCACCTGG-GCCCCTAAGCCAGA-TCTAGGGAATCCCT-GCATGGCAAACCGCCACTAACCC-TCAGTTTTTGTCAAGCTGTCAGA-CCTCCAGCGCAGGTTTCAGTGCCATTCATGTC-CCTGCGAGTGCTTATCAATGGTT-TATGACGGATATCCCACATTCGGAGAACA-AAACAGGAGAA-GATCTTGAATATGGGGCATGTCCTAATAACATGATGGGCAC-TTCTCAGTGCGGACTGTGGGGACCTCCAAGTCCAAGTACCCTTTAGTGGTTAGGATTTACATGAGAATGAAGCACGTCAGGGCGTGGATACCTCGCCCGATGCG-AACCAGAA-TACCTATTCAAAGCCAACCCAAATTATGC-GGCAACTCCATTAAGCCAACTGGT-CCAGTCGCACAGC-ATCACTACTCTT";
            String tmpString2 = "GAGATAGGGTGGCAGATGTAATTGAAAGTTCCATAGGAGATAGCGTG-GCAGAGCCCTCACTCACGCTCTACCAGCACCCACAGGCCAG-ACACACAGGTGAGCAGTCATCGACTGGA-ACAGGCAAGGTTCCAGCACTCCAAGCTGCTGAAATTGGAGCATCATCAAATGCTAGTGACGAGAG-ATGATTGAGACACGCTGTGT-CTTAACTCGCACAGCACAGCTGAGAC-ACTCTTGATAGTTTCTTCAGCAGAGCGGGATTAGTTGGAGAG-TAGATCTCCCTCTTGAAGGCACAACTAACCCAAATGGTTATGCCAACTGGGACATAGATATAACAGGTTACGCGCAAATGCGTAGAAAGGTAGAGCT-TTCACCTACATGCGCTTTGATGCAGAGTTCACTTTTGTTGCGTGCACACCCACCGGG--AGTTGTCCCACAATTGCTCCAATATATGTTTGTGCCACCTGG-GCCCCTAAGCCAGA-TCTAGGGAATCCCT-GCATGGCAAACCGCCACTAACCC-TCAGTTTTTGTCAAGCTGTCAGA-CCTCCAGCGCAGGTTTCAGTGCCATTCATGTC-CCTGCGAGTGCTTATCAATGGTT-TATGACGGATATCCCACATTCGGAGAACA-AAACAGGAGAA-GATCTTGAATATGGGGCATGTCCTAATAACATGATGGGCAC-TTCTCAGTGCGGACTGTGGGGACCTCCAAGTCCAAGTACCCTTTAGTGGTTAGGATTTACATGAGAATGAAGCACGTCAGGGCGTGGATACCTCGCCCGATGCG-AACCAGAA-TACCTATTCAAAGCCAACCCAAATTATGC-GGCAACTCCATTAAGCCAACTGGT-CCAGTCGCACAGC-ATCACTACTCTT";
            
            String tmpString3 = "GAGATAGGGTGGCAGATGTAATTGAAAGTTCCATAGGAGATAGCGTG-GCAGAGCCCTCACTCACGCTCTACCAGCACCCACAGGCCAG-ACACACAGGTGAGCAGTCATCGACTGGA-ACAGGCAAGGTTCCAGCACTCCAAGCTGCTGAAATTGGAGCATCATCAAATGCTAGTGACGAGAG-ATGATTGAGACACGCTGTGT-CTTAACTCGCACAGCACAGCTGAGAC-ACTCTTGATAGTTTCTTCAGCAGAGCGGGATTAGTTGGAGAG-TAGATCTCCCTCTTGAAGGCACAACTAACCCAAATGGTTATGCCAACTGGGACATAGATATAACAGGTTACGCGCAAATGCGTAGAAAGGTAGAGCT-TTCACCTACATGCGCTTTGATGCAGAGTTCACTTTTGTTGCGTGCACACCCACCGGG--AGTTGTCCCACAATTGCTCCAATATATGTTTGTGCCACCTGG-GCCCCTAAGCCAGA-TCTAGGGAATCCCT-GCATGGCAAACCGCCACTAACCC-TCAGTTTTTGTCAAGCTGTCAGA-CCTCCAGCGCAGGTTTCAGTGCCATTCATGTC-CCTGCGAGTGCTTATCAATGGTT-TATGACGGATATCCCACATTCGGAGAACA-AAACAGGAGAA-GATCTTGAATATGGGGCATGTCCTAATAACATGATGGGCAC-TTCTCAGTGCGGACTGTGGGGACCTCCAAGTCCAAGTACCCTTTAGTGGTTAGGATTTACATGAGAATGAAGCACGTCAGGGCGTGGATACCTCGCCCGATGCG-AACCAGAA-TACCTATTCAAAGCCAACCCAAATTATGC-GGCAACTCCATTAAGCCAACTGGT-CCAGTCGCACAGC-ATCACTACTCTT";
            nrSeqs[0] = DNATools.createGappedDNASequence(tmpString1, "Seq1");
            nrSeqs[1] = DNATools.createGappedDNASequence(tmpString2, "Seq2");
            nrSeqs[2] = DNATools.createGappedDNASequence(tmpString3, "Seq3");
            

            aaSeqs[0] = new ProteinSequence("IFRGLEICCYGPFTNMPTDQLEWMVQLCGASVVKELSSFTLGTGVHPIVVVQPDAWT---GFHAIGQMCEAPVVTREWVLDSVALYQCQELDTYLIPQIP-------");
            AccessionID accId4 = new AccessionID("AA22551");
            aaSeqs[0].setAccession(accId4);

            aaSeqs[1] = new ProteinSequence("IFRGLEICCYGPFTNMPTDQLEWMVQLCGASVVKELSSFTLGTGVHPIVVVQPDAWT---GFHAIGQMCEAPVVTREWVLDSVALYQCQELDTYLIPQIP-------");
            AccessionID accId5 = new AccessionID("BB22551");
            aaSeqs[1].setAccession(accId5);

            aaSeqs[2] = new ProteinSequence("IFRGLEICCYGPFTNMPTDQLEWMVQLCGASVVKELSSFTLGTGVHPIVVVQPDAWT---GFHAIGQMCEAPVVTREWVLDSVALYQCQELDTYLIPQIP-------");
            AccessionID accId6 = new AccessionID("CC22551");
            aaSeqs[2].setAccession(accId6);
        } 
        catch (BioException bioe) 
        {
            System.err.println("Bioexception: " + bioe);
        }
        
        SeqPaint seqPaint = new SeqPaint();
        seqPaint.paintSeq(nrSeqs);
//        seqPaint.paintSeq(aaSeqs);

    }

    private static void testParsStructure() throws FileNotFoundException, IOException 
    {
        String filePath = "/home/zhongxf/atest/test/data0407.pdb";
        String filePath1 = "/home/zhongxf/aTEST/BRCA1/1JNX.pdb";
        String filePath2 = "/home/zhongxf/aTEST/BRCA1/1N5O.pdb";
        File file = new File(filePath);
        File[] files = {new File(filePath1), new File(filePath2) };
        Structure[] structArray = new Structure[files.length];
         for (int i =0; i< files.length; i ++)
         {
             String fileName = files[i].getName();
            FileReader fileReader = new FileReader(files[i]);
            BufferedReader buffReader = new BufferedReader(fileReader);
            PDBFileReader pdbRead = new PDBFileReader();
            structArray[i] = pdbRead.getStructure(files[i]); 
         }
         
        ParseStructure parsSt = new ParseStructure();
        //assume the structArray is a collection of structure alignment
        parsSt.ParsStruct(structArray);
        parsSt.ParsStruct(parsSt.parsFile(file));
        //extract structures from an alignment file generated by structure alignment tools;and return structure[]
        parsSt.parsFile(file);
    }

   
}
