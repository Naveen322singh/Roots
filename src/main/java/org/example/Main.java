package org.example;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.JsonNode;
import java.io.File;
import java.math.BigInteger;
import java.util.*;

public class Main {

    // A simple container for one root (x, y)
    static class Root {
        double x;              // x-coordinate
        double y;              // y-coordinate (converted from given base)
        int base;              // base of the input number
        String originalValue;  // original number as string

        Root(double x, double y, int base, String originalValue) {
            this.x = x;
            this.y = y;
            this.base = base;
            this.originalValue = originalValue;
        }
    }

    // Convert a number given in some base into decimal
    static double convertToDecimal(String number, int base) {
        try {
            return new BigInteger(number.toLowerCase(), base).doubleValue();
        } catch (Exception e) {
            return 0; // if something goes wrong, return 0
        }
    }

    // Read the JSON file and turn it into a list of roots
    static List<Root> readRoots(String filename) throws Exception {
        ObjectMapper mapper = new ObjectMapper();
        JsonNode rootNode = mapper.readTree(new File(filename));

        List<Root> roots = new ArrayList<>();

        // go through each entry in the JSON
        Iterator<Map.Entry<String, JsonNode>> fields = rootNode.fields();
        while (fields.hasNext()) {
            Map.Entry<String, JsonNode> entry = fields.next();
            String key = entry.getKey();

            // skip the "keys" object
            if ("keys".equals(key)) continue;

            try {
                double x = Double.parseDouble(key); // the key itself is x
                JsonNode rootData = entry.getValue();
                int base = rootData.get("base").asInt();
                String valueStr = rootData.get("value").asText();
                double y = convertToDecimal(valueStr, base);

                roots.add(new Root(x, y, base, valueStr));
            } catch (NumberFormatException ignored) {}
        }

        // sort them by x value
        roots.sort(Comparator.comparingDouble(r -> r.x));
        return roots;
    }

    // Read k (how many points are needed) from JSON
    static int readK(String filename) throws Exception {
        ObjectMapper mapper = new ObjectMapper();
        JsonNode rootNode = mapper.readTree(new File(filename));
        return rootNode.get("keys").get("k").asInt();
    }

    // Build the Vandermonde system and solve it to get polynomial coefficients
    static double[] findPolynomialCoefficients(List<Root> points, int degree) {
        int n = points.size();
        double[][] matrix = new double[n][degree + 2];

        for (int i = 0; i < n; i++) {
            Root r = points.get(i);
            for (int j = 0; j <= degree; j++) {
                matrix[i][j] = Math.pow(r.x, j);
            }
            matrix[i][degree + 1] = r.y; // RHS (y value)
        }

        return solveByGaussianElimination(matrix);
    }

    // Solve a system of equations with Gaussian elimination
    static double[] solveByGaussianElimination(double[][] matrix) {
        int n = matrix.length;
        int m = matrix[0].length - 1;

        // forward elimination
        for (int i = 0; i < n; i++) {
            // pivot: find the row with the biggest element in this column
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(matrix[k][i]) > Math.abs(matrix[maxRow][i])) {
                    maxRow = k;
                }
            }
            // swap rows
            double[] temp = matrix[i];
            matrix[i] = matrix[maxRow];
            matrix[maxRow] = temp;

            // eliminate below
            for (int k = i + 1; k < n; k++) {
                double factor = matrix[k][i] / matrix[i][i];
                for (int j = i; j <= m; j++) {
                    matrix[k][j] -= factor * matrix[i][j];
                }
            }
        }

        // back substitution
        double[] solution = new double[m];
        for (int i = n - 1; i >= 0; i--) {
            solution[i] = matrix[i][m];
            for (int j = i + 1; j < m; j++) {
                solution[i] -= matrix[i][j] * solution[j];
            }
            solution[i] /= matrix[i][i];
        }

        return solution;
    }

    public static void main(String[] args) {
        try {
            String filename = "testcase.json";

            // read input
            List<Root> roots = readRoots(filename);
            int k = readK(filename);

            // not enough points to build polynomial
            if (roots.size() < k) {
                System.out.println("Error: Not enough roots to solve");
                return;
            }

            // use only the first k points
            List<Root> selected = roots.subList(0, k);

            // solve for coefficients
            double[] coefficients = findPolynomialCoefficients(selected, k - 1);

            // constant term (c) is the first coefficient
            double c = coefficients[0];

            // print just the constant
            System.out.printf("%.0f%n", c);

        } catch (Exception e) {
            System.err.println("Something went wrong: " + e.getMessage());
        }
    }
}
