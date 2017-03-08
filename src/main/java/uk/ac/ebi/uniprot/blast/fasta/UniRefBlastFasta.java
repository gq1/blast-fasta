package uk.ac.ebi.uniprot.blast.fasta;

import uk.ac.ebi.uniprot.blast.fasta.FastaReader.Sequence;
import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.*;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.input.DatabaseOption;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class provides a simple example of submitting blast jobs through the UniProtJAPI.
 */
public class UniRefBlastFasta {

    private static final Logger logger = LoggerFactory.getLogger(UniRefBlastFasta.class);
    private static int numberOfThread = 1;

    public static void main(String[] args) throws IOException {
        String fastaFilename = args[0];
        runUniRefBlast(fastaFilename, DatabaseOption.UNIREF_90);
    }

    public static void runUniRefBlast(String fastaFile, DatabaseOption uniRefDB) throws IOException {

        logger.info("Start UniRef blast");

        ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
        UniRefBlastService uniRefBlastService = serviceFactoryInstance.getUniRefBlastService();
        uniRefBlastService.start();

        LimitedQueue<Runnable> queue = new LimitedQueue<>(10000);

        int coreNumber = Runtime.getRuntime().availableProcessors();
        ExecutorService executorService =
                new ThreadPoolExecutor(Math.max(numberOfThread, coreNumber), 32, 60, TimeUnit.SECONDS, queue);
        
        PrintWriter outputWriter
        = new PrintWriter(new BufferedWriter(new FileWriter("fastaBlast.tsv")));
        
        FastaReader fastaReader = new FastaReader(fastaFile);
        fastaReader.openInputFile();
        FastaReader.Sequence sequence;
        while ((sequence = fastaReader.nextSequence()) != null) {
            Runnable runnable = new BlastRunnable(sequence, uniRefDB, uniRefBlastService, outputWriter);
            executorService.submit(runnable);
        }

        executorService.shutdown();
        logger.info("Shutting down executor");

        try {
            executorService.awaitTermination(15, TimeUnit.MINUTES);
        } catch (InterruptedException e) {
            logger.warn("Executer waiting interrupted.", e);
        }

        fastaReader.close();
        uniRefBlastService.stop();
        outputWriter.close();
        logger.info("Finished UniRef blast");
    }

    static public class LimitedQueue<E> extends LinkedBlockingQueue<E> {

        private static final long serialVersionUID = -1886696548851000684L;

        public LimitedQueue(int maxSize) {
            super(maxSize);
        }

        @Override
        public boolean offer(E e) {
            // turn offer() and add() into a blocking calls (unless interrupted)
            try {
                put(e);
                return true;
            } catch (InterruptedException ie) {
                Thread.currentThread().interrupt();
            }
            return false;
        }
    }

    public static class BlastRunnable implements Runnable {

        private static int count = 0;
        private Sequence sequence;
        private DatabaseOption uniRefDB;
        private UniRefBlastService uniRefBlastService;
        private PrintWriter outputWriter;

        public BlastRunnable(Sequence sequence, DatabaseOption uniRefDB, UniRefBlastService uniRefBlastService, PrintWriter outputWriter) {

            this.sequence = sequence;
            this.uniRefDB = uniRefDB;
            this.uniRefBlastService = uniRefBlastService;
            this.outputWriter = outputWriter;
        }

        @Override
        public void run() {
            count++;
            logger.info("Blasting {}: {}", count, sequence.name);

            BlastInput input = new BlastInput.Builder(uniRefDB, sequence.seq).build();

            CompletableFuture<BlastResult<UniRefHit>> resultFuture = uniRefBlastService.runBlast(input);

            try {
                BlastResult<UniRefHit> blastResult = resultFuture.get();
                logger.info("Number of blast hits: " + blastResult.getNumberOfHits());

                for (UniRefHit hit : blastResult.hits()) {
               
                   outputWriter.println(
                           sequence.name + "\t"
                         + hit.getEntry().getUniRefEntryId().getValue() + "\t"
                         + hit.getSummary()                      
                           );                   
                    break;
                }
            } catch (ExecutionException e) {
                logger.error(e.getCause().getMessage());
            } catch (InterruptedException e) {
                logger.error(e.getMessage());
            }
        }

    }
}