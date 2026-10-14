import java.util.*;
class MatrixException extends Exception {
    MatrixException(String s) {
        super(s);
    }
}

class Complex {
    private double re, im;

    Complex() {
        re = 0;
        im = 0;
    }

    Complex(double re, double im) {
        this.re = re;
        this.im = im;
    }

    public double getRe() {
        return re;
    }

    public double getIm() {
        return im;
    }

    public void setRe(double re) {
        this.re = re;
    }

    public void setIm(double im) {
        this.im = im;
    }

    Complex slozh(Complex b) {
        return new Complex(re + b.re, im + b.im);
    }

    Complex vych(Complex b) {
        return new Complex(re - b.re, im - b.im);
    }

    Complex umnozh(Complex b) {
        return new Complex(re * b.re - im * b.im, re * b.im + im * b.re);
    }

    Complex delit(Complex b) {
        double z = b.re * b.re + b.im * b.im;
        if (z == 0) {
            throw new ArithmeticException("деление на ноль");
        }
        return new Complex((re * b.re + im * b.im) / z, (im * b.re - re * b.im) / z);
    }

    boolean nol() {
        return re == 0 && im == 0;
    }

    public String toString() {
        double a = Math.round(re * 1000) / 1000.0;
        double b = Math.round(im * 1000) / 1000.0;
        if (b < 0) {
            return a + "-" + (-b) + "i";
        }
        return a + "+" + b + "i";
    }
}

class Matrix {
    private int n, k;
    private Complex m[][];

    Matrix(int n, int k) {
        this.n = n;
        this.k = k;
        m = new Complex[n][k];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < k; j++) {
                m[i][j] = new Complex();
            }
        }
    }

    public int getN() {
        return n;
    }

    public int getK() {
        return k;
    }

    public Complex get(int i, int j) {
        return m[i][j];
    }

    public void set(int i, int j, Complex c) {
        m[i][j] = c;
    }

    Matrix slozh(Matrix b) throws MatrixException {
        if (n != b.n || k != b.k) {
            throw new MatrixException("размеры матриц не совпадают");
        }
        Matrix rez = new Matrix(n, k);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < k; j++) {
                rez.m[i][j] = m[i][j].slozh(b.m[i][j]);
            }
        }
        return rez;
    }

    Matrix umnozh(Matrix b) throws MatrixException {
        if (k != b.n) {
            throw new MatrixException("число столбцов первой не равно числу строк второй");
        }
        Matrix rez = new Matrix(n, b.k);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < b.k; j++) {
                Complex sum = new Complex();
                for (int t = 0; t < k; t++) {
                    sum = sum.slozh(m[i][t].umnozh(b.m[t][j]));
                }
                rez.m[i][j] = sum;
            }
        }
        return rez;
    }

    Matrix transp() {
        Matrix rez = new Matrix(k, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < k; j++) {
                rez.m[j][i] = m[i][j];
            }
        }
        return rez;
    }

    Matrix minor(int a, int b) {
        Matrix rez = new Matrix(n - 1, k - 1);
        int x = 0;
        for (int i = 0; i < n; i++) {
            if (i == a) continue;
            int y = 0;
            for (int j = 0; j < k; j++) {
                if (j == b) continue;
                rez.m[x][y] = m[i][j];
                y++;
            }
            x++;
        }
        return rez;
    }

    Complex opred() throws MatrixException {
        if (n != k) {
            throw new MatrixException("матрица не квадратная");
        }
        if (n == 1) {
            return m[0][0];
        }
        if (n == 2) {
            return m[0][0].umnozh(m[1][1]).vych(m[0][1].umnozh(m[1][0]));
        }
        Complex rez = new Complex();
        int znak = 1;
        for (int j = 0; j < k; j++) {
            Complex c = m[0][j].umnozh(minor(0, j).opred());
            if (znak == 1) {
                rez = rez.slozh(c);
            } else {
                rez = rez.vych(c);
            }
            znak = -znak;
        }
        return rez;
    }

    Matrix obratn() throws MatrixException {
        Complex d = opred();
        if (d.nol()) {
            throw new MatrixException("определитель равен нулю, обратной нет");
        }
        Matrix rez = new Matrix(n, n);
        if (n == 1) {
            rez.m[0][0] = new Complex(1, 0).delit(d);
            return rez;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Complex c = minor(i, j).opred();
                if ((i + j) % 2 == 1) {
                    c = new Complex().vych(c);
                }
                rez.m[j][i] = c.delit(d);
            }
        }
        return rez;
    }

    Matrix delit(Matrix b) throws MatrixException {
        if (b.n != b.k) {
            throw new MatrixException("делить можно только на квадратную матрицу");
        }
        return umnozh(b.obratn());
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append("[ ");
            for (int j = 0; j < k; j++) {
                sb.append(m[i][j] + " ");
            }
            sb.append("]\n");
        }
        return sb.toString();
    }
}

public class Main {
    static Matrix vvod(Scanner sc) throws MatrixException {
        System.out.println("строки:");
        int n = sc.nextInt();
        System.out.println("столбцы:");
        int k = sc.nextInt();
        if (n <= 0 || k <= 0) {
            throw new MatrixException("неправильный размер");
        }
        Matrix rez = new Matrix(n, k);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < k; j++) {
                System.out.println("[" + i + "][" + j + "] re im:");
                double x = sc.nextDouble();
                double y = sc.nextDouble();
                rez.set(i, j, new Complex(x, y));
            }
        }
        return rez;
    }

    public static void main(String args[]) {
        Scanner sc = new Scanner(System.in);
        Matrix a = null;
        Matrix b = null;
        int vybor = -1;
        while (vybor != 0) {
            System.out.println("1 - ввод a");
            System.out.println("2 - ввод b");
            System.out.println("3 - вывод a и b");
            System.out.println("4 - a+b");
            System.out.println("5 - a*b");
            System.out.println("6 - a/b");
            System.out.println("7 - транспонировать a");
            System.out.println("8 - транспонировать b");
            System.out.println("9 - определитель a");
            System.out.println("10 - определитель b");
            System.out.println("0 - выход");
            vybor = sc.nextInt();
            try {
                switch (vybor) {
                    case 1:
                        a = vvod(sc);
                        break;
                    case 2:
                        b = vvod(sc);
                        break;
                    case 3:
                        System.out.println("a:");
                        System.out.println(a);
                        System.out.println("b:");
                        System.out.println(b);
                        break;
                    case 4:
                        if (a == null || b == null) throw new MatrixException("нет матрицы");
                        System.out.println(a.slozh(b));
                        break;
                    case 5:
                        if (a == null || b == null) throw new MatrixException("нет матрицы");
                        System.out.println(a.umnozh(b));
                        break;
                    case 6:
                        if (a == null || b == null) throw new MatrixException("нет матрицы");
                        System.out.println(a.delit(b));
                        break;
                    case 7:
                        if (a == null) throw new MatrixException("нет матрицы");
                        a = a.transp();
                        System.out.println(a);
                        break;
                    case 8:
                        if (b == null) throw new MatrixException("нет матрицы");
                        b = b.transp();
                        System.out.println(b);
                        break;
                    case 9:
                        if (a == null) throw new MatrixException("нет матрицы");
                        System.out.println("det = " + a.opred());
                        break;
                    case 10:
                        if (b == null) throw new MatrixException("нет матрицы");
                        System.out.println("det = " + b.opred());
                        break;
                    case 0:
                        break;
                    default:
                        System.out.println("нет такого");
                }
            }
            catch (MatrixException e) {
                System.out.println("ошибка: " + e.getMessage());
            }
            catch (ArithmeticException e) {
                System.out.println("ошибка: " + e.getMessage());
            }
        }
    }
}
