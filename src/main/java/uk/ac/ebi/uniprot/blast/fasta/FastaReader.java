package uk.ac.ebi.uniprot.blast.fasta;


import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Map;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author gqi
 *
 *         This class is to read a fasta file sequence by sequence
 *
 */
public class FastaReader implements Closeable {

    private static final Logger log = LoggerFactory.getLogger(FastaReader.class);

    private final File file;
    private BufferedReader fileReader;
    private String currentLine;
    private int currentNumberOfSequences;

    /**
     *
     * @param fileName
     * @throws IOException
     * @throws FileNotFoundException
     */
    public FastaReader(String fileName) {

        this.file = new File(fileName);
    }

    /**
     *
     * @param file
     * @throws FileNotFoundException
     * @throws IOException
     */

    public FastaReader(File file) {
        this.file = file;
    }

    @Override
    public void close() {

        if (this.fileReader != null) {
            try {
                fileReader.close();
            } catch (IOException e) {
                log.error("Cannot close fastq file reader", e);
            }
        }
    }

    /**
     *
     * @return next sequence in fasta file or null when reach the end of file
     */
    public Sequence nextSequence() {

        if (this.currentLine == null) {

            log.debug("There is no more sequence in the fasta file");
            return null;
        } else if (!this.currentLine.startsWith(">")) {

            log.error("The forma of fasta file is wrong");
        } else {

            String name = this.currentLine.substring(1);
            StringBuilder seqBuilder = new StringBuilder();

            do {

                try {
                    this.currentLine = this.fileReader.readLine();
                    if (currentLine != null && !currentLine.startsWith(">")) {
                        seqBuilder.append(currentLine);
                    }
                } catch (IOException e) {
                    log.error("Problems to read fasta file", e);
                }
            } while (this.currentLine != null
                    && !this.currentLine.startsWith(">"));

            this.currentNumberOfSequences++;

            String seq = seqBuilder.toString();
            return new Sequence(name, seq);
        }
        return null;
    }

    /**
     *
     * @param file
     * @throws FileNotFoundException
     * @throws IOException
     */
    public void openInputFile() throws IOException {

        String fileName = file.getPath();

        if (!file.exists()) {

            throw new FileNotFoundException("File does not exist: " + fileName);
        } else if (file.isDirectory()) {

            throw new IllegalArgumentException("File name is a directory: "
                    + fileName);
        } else if (!file.canRead()) {

            throw new FileNotFoundException("File cannot be read: " + fileName);
        } else {

            this.fileReader = new BufferedReader(new FileReader(file));
            this.currentLine = fileReader.readLine();
        }

    }

    /**
     * @return the currentNumberOfSequences
     */
    public int getCurrentNumberOfSequences() {

        return currentNumberOfSequences;
    }

    /**
     *
     * @return
     */
    public BufferedReader getFileReader() {
        return fileReader;
    }

    /**
     *
     * @author gqi
     *
     *         This class represents a sequence with name and sequence
     *
     */
    public static class Sequence {

        public String name;
        public String seq;

        public Sequence(String name, String seq) {
            this.name = name;
            this.seq = seq;
        }

        /**
         *
         * @return
         */
        public int getSeqLength() {
            return seq.length();
        }

        /**
         *
         * @return
         */
        public String getMd5() {

            StringBuilder result = new StringBuilder();

            MessageDigest complete;
            try {
                complete = MessageDigest.getInstance("MD5");
            } catch (NoSuchAlgorithmException e) {
                log.error("Problems to calculate md5", e);
                return null;
            }

            complete.update(this.seq.getBytes());
            byte[] b = complete.digest();

            for (int i = 0; i < b.length; i++) {
                result.append(Integer.toString((b[i] & 0xff) + 0x100, 16)
                        .substring(1));
            }
            return result.toString();
        }

        /**
         *
         * @return
         */
        public Map<Character, Integer> getAAComposition() {

            HashMap<Character, Integer> aaComposition = new HashMap<>();

            for (char aa : this.seq.toUpperCase().toCharArray()) {
                if (aa == '*') {
                    continue;
                }
                Integer freq = aaComposition.get(aa);
                if (freq == null) {
                    aaComposition.put(aa, 1);
                } else {
                    aaComposition.put(aa, freq.intValue() + 1);
                }
            }

            return aaComposition;
        }

    }

}
