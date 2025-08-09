import java.io.*;
import java.math.BigInteger;
import java.util.*;
import org.json.JSONObject; 

public class PolynomialConstantFinder {

    // Helper class for exact rationals
    static class BigFraction {
        BigInteger num, den;

        BigFraction(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new ArithmeticException("Divide by zero");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            num = n.divide(g);
            den = d.divide(g);
        }

        static BigFraction fromInt(BigInteger n) { return new BigFraction(n, BigInteger.ONE); }

        BigFraction add(BigFraction o) {
            return new BigFraction(
                    num.multiply(o.den).add(o.num.multiply(den)),
                    den.multiply(o.den)
            );
        }

        BigFraction subtract(BigFraction o) {
            return new BigFraction(
                    num.multiply(o.den).subtract(o.num.multiply(den)),
                    den.multiply(o.den)
            );
        }

        BigFraction multiply(BigFraction o) {
            return new BigFraction(num.multiply(o.num), den.multiply(o.den));
        }

        BigFraction divide(BigFraction o) {
            return new BigFraction(num.multiply(o.den), den.multiply(o.num));
        }
    }

    // Decode string value in given base to BigInteger
    static BigInteger decode(String value, int base) {
        BigInteger result = BigInteger.ZERO;
        BigInteger b = BigInteger.valueOf(base);
        for (char ch : value.toCharArray()) {
            int digit;
            if (ch >= '0' && ch <= '9') digit = ch - '0';
            else if (ch >= 'a' && ch <= 'z') digit = 10 + (ch - 'a');
            else if (ch >= 'A' && ch <= 'Z') digit = 10 + (ch - 'A');
            else throw new IllegalArgumentException("Invalid digit: " + ch);
            result = result.multiply(b).add(BigInteger.valueOf(digit));
        }
        return result;
    }

    // Solve linear system Ax = y for polynomial coefficients
    static BigFraction[] solveVandermonde(BigInteger[] xs, BigInteger[] ys) {
        int n = xs.length;
        BigFraction[][] A = new BigFraction[n][n];
        BigFraction[] Y = new BigFraction[n];

        // Build Vandermonde matrix
        for (int i = 0; i < n; i++) {
            BigInteger pow = BigInteger.ONE;
            for (int j = n - 1; j >= 0; j--) {
                A[i][n - 1 - j] = BigFraction.fromInt(xs[i].pow(j));
            }
            Y[i] = BigFraction.fromInt(ys[i]);
        }

        // Gaussian elimination
        for (int col = 0; col < n; col++) {
            // Find pivot
            int pivot = col;
            for (int i = col; i < n; i++) {
                if (!A[i][col].num.equals(BigInteger.ZERO)) { pivot = i; break; }
            }
            if (pivot != col) {
                BigFraction[] tmp = A[col]; A[col] = A[pivot]; A[pivot] = tmp;
                BigFraction t = Y[col]; Y[col] = Y[pivot]; Y[pivot] = t;
            }

            // Normalize pivot row
            BigFraction inv = new BigFraction(A[col][col].den, A[col][col].num); // reciprocal
            for (int j = col; j < n; j++) A[col][j] = A[col][j].multiply(inv);
            Y[col] = Y[col].multiply(inv);

            // Eliminate below
            for (int i = col + 1; i < n; i++) {
                BigFraction factor = A[i][col];
                if (factor.num.equals(BigInteger.ZERO)) continue;
                for (int j = col; j < n; j++)
                    A[i][j] = A[i][j].subtract(factor.multiply(A[col][j]));
                Y[i] = Y[i].subtract(factor.multiply(Y[col]));
            }
        }

        // Back substitution
        BigFraction[] coeff = new BigFraction[n];
        for (int i = n - 1; i >= 0; i--) {
            BigFraction sum = BigFraction.fromInt(BigInteger.ZERO);
            for (int j = i + 1; j < n; j++)
                sum = sum.add(A[i][j].multiply(coeff[j]));
            coeff[i] = Y[i].subtract(sum);
        }
        return coeff;
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: java PolynomialConstantFinder <input.json>");
            return;
        }
        String jsonText = new String(java.nio.file.Files.readAllBytes(java.nio.file.Paths.get(args[0])));
        JSONObject obj = new JSONObject(jsonText);

        int n = obj.getJSONObject("keys").getInt("n");
        int k = obj.getJSONObject("keys").getInt("k");

        List<Integer> xs = new ArrayList<>();
        List<BigInteger> ys = new ArrayList<>();

        // Read & decode all roots
        for (String key : obj.keySet()) {
            if (key.equals("keys")) continue;
            int x = Integer.parseInt(key);
            JSONObject root = obj.getJSONObject(key);
            int base = Integer.parseInt(root.getString("base"));
            String val = root.getString("value");
            BigInteger y = decode(val, base);
            xs.add(x);
            ys.add(y);
        }

        // Sort by x and pick first k roots
        List<Integer> sortedIdx = new ArrayList<>();
        for (int i = 0; i < xs.size(); i++) sortedIdx.add(i);
        sortedIdx.sort(Comparator.comparingInt(xs::get));

        BigInteger[] selX = new BigInteger[k];
        BigInteger[] selY = new BigInteger[k];
        for (int i = 0; i < k; i++) {
            selX[i] = BigInteger.valueOf(xs.get(sortedIdx.get(i)));
            selY[i] = ys.get(sortedIdx.get(i));
        }

        // Solve
        BigFraction[] coeffs = solveVandermonde(selX, selY);

        // Constant term is last coefficient (a0)
        BigFraction cFrac = coeffs[coeffs.length - 1];
        if (!cFrac.den.equals(BigInteger.ONE)) {
            System.out.println("c is not an integer: " + cFrac.num + "/" + cFrac.den);
        } else {
            System.out.println("c = " + cFrac.num);
        }
    }
}
