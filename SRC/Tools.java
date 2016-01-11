package ARGON;

import java.util.Random;

/**
 *
 * @author Pier Palamara <ppalama@hsph.harvard.edu>
 */
public class Tools {

    public static double sampleExponential(double rate, Random generator) {
        return -Math.log(generator.nextDouble()) / rate;
    }

    static int roundToBlock(int length, int blockSize) {
        int excess = length % blockSize;
        return length - excess + (int) Math.round((double) excess / blockSize) * blockSize;
    }

// using Lanczos approximation of Gamma. See e.g.
//    https://en.wikipedia.org/wiki/Lanczos_approximation
//    http://mrob.com/pub/ries/lanczos-gamma.html
//    Java translation from c++ implementation by Liming Liang and Goncalo Abecasis (GENOME: http://csg.sph.umich.edu/liang/genome/)
    
 static int samplePoisson(double lambda, Random generator) {

        double waitTime, s, logmean, c;
        double count, eventTime, ratio, y;
        double maxInt = Integer.MAX_VALUE;
        if (lambda < 12) {
            waitTime = Math.exp(-lambda);

            count = 0;
            eventTime = generator.nextDouble();

            while (eventTime > waitTime) {
                count++;
                eventTime *= generator.nextDouble();
            }
        } else {
            s = Math.sqrt(2 * lambda);
            logmean = Math.log(lambda);
            c = lambda * Math.log(lambda) - lngamma(lambda + 1);
            do {
                do {
                    y = Math.tan(Math.PI * generator.nextDouble());
                    count = Math.floor(s * y + lambda);
                } while (count < 0);

                ratio = 0.9 * (1 + y * y) * Math.exp(count * logmean - lngamma(count + 1) - c);

            } while (generator.nextDouble() > ratio);
        }
        if (count > maxInt) {
            count = maxInt;
        }
        return (int) count;
        
    }

    static double lngamma(double z) {
        double result, sum;
        double[] c = new double[7];
        c[0] = 1.000000000190015;
        c[1] = 76.18009172947146;
        c[2] = -86.50532032941677;
        c[3] = 24.01409824083091;
        c[4] = -1.231739572450155;
        c[5] = 0.1208650973866179e-2;
        c[6] = -0.5395239384953e-5;
        sum = c[0];
        for (int i = 1; i < c.length; i++) {
            sum += c[i] / (z + i);
        }
        result = (z + 0.5) * Math.log(z + 5.5) - (z + 5.5) + Math.log(2.5066282746310005 * sum / z);
        return result;
    }

}
