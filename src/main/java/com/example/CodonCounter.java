package com.example;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.Map;

public class CodonCounter {

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: CodonCounter <genome.fasta>");
            System.exit(1);
        }

        LinkedHashMap<String, DNASequence> records =
                FastaReaderHelper.readFastaDNASequence(new File(args[0]));

        Map<String, Integer> counts = new LinkedHashMap<>();
        long totalWindows = 0;

        for (DNASequence record : records.values()) {
            String seq = record.getSequenceAsString().toUpperCase();
            for (int i = 0; i + 3 <= seq.length(); i++) {
                String codon = seq.substring(i, i + 3);
                if (!codon.matches("[ACGT]{3}")) continue;
                counts.merge(codon, 1, Integer::sum);
                totalWindows++;
            }
        }

        TranscriptionEngine engine = new TranscriptionEngine.Builder()
                .initMet(false)
                .trimStop(false)
                .build();

        System.out.printf("Genome records read: %d | Windows counted: %d%n%n", records.size(), totalWindows);
        System.out.printf("%-4s %-7s %-7s %-10s %-9s %s%n",
                "#", "DNA", "RNA", "Count", "Freq %", "Amino acid");

        final long total = totalWindows;
        int[] rank = {1};
        counts.entrySet().stream()
                .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
                .limit(10)
                .forEach(e -> {
                    DNASequence dna = null;
                    try {
                        dna = new DNASequence(e.getKey());
                        RNASequence rna = dna.getRNASequence(engine);
                        ProteinSequence aa = rna.getProteinSequence(engine);
                        String letter = aa.getSequenceAsString();
                        if (letter.isEmpty() || letter.equals("*")) letter = "Stop";
                        System.out.printf("%-4d %-7s %-7s %-10d %-9.3f %s%n",
                                rank[0]++, e.getKey(), rna.getSequenceAsString(),
                                e.getValue(), 100.0 * e.getValue() / total, letter);
                    } catch (Exception ex) {
                        throw new RuntimeException(ex);
                    }
                });
    }
}