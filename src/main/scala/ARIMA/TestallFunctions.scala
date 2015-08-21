/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

// package TimeSeriesForecast.ARIMA

package ARIMA

import breeze.linalg._
import breeze.stats._
import breeze.numerics._
import breeze.interpolation.LinearInterpolator
import org.apache.commons.math3.distribution.ChiSquaredDistribution
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression
import stat.StatDenseVector
import IO.Readcsv._

// this is an object to test most of functions used in the implementations
// should not be used
object TestallFunctions {

  def toeplitz(x: DenseVector[Double]) = {
    val l = x.length
    var A = DenseMatrix.zeros[Double](l, l)
    A = new DenseMatrix[Double](l, l, (abs(col(A) - row(A))).map(v => x(v)).toArray)
    A
  }

  def ChampernowneD(z: DenseVector[Double], p: Int, MeanZero: Boolean = false): DenseMatrix[Double] = {
    val n = z.length
    var y = DenseVector.zeros[Double](n)
    if (MeanZero) y = z
    else y = z - mean(z)
    var x = y.toArray
    val x0 = y
    for (i <- 0 until p) {
      x = x ++ (Array.fill(i + 1)(0.0) ++ x.slice(0, n - i - 1))
    }
    val newx = DenseMatrix.zeros[Double](p + 1, x0.length)
    for (i <- 0 to p) {
      for (j <- 0 until x0.length) {
        newx(i, j) = x(i * x0.length + j)
      }
    }
    val C = newx * x0
    val A = toeplitz(C)
    val E = DenseMatrix.zeros[Double](p + 1, p + 1)
    for (j <- 0 until p) for (i <- 0 to j) E(i + 1, j + 1) = E(i, j) + y(i) * y(j) + y(n - 1 - i) * y(n - 1 - j)
    for (j <- 0 to p)
      for (i <- 0 to j - 1)
        E(j, i) = E(i, j)
    A - E
  }

  def PacfToAR(zeta: DenseVector[Double]): DenseVector[Double] = {
    val L = zeta.length
    if (L == 0) return (DenseVector.zeros[Double](0))
    if (L == 1) return (zeta)
    var phik = new DenseVector[Double](Array(zeta(0)))
    for (k <- 1 until L) {
      var phikm1 = phik
      phik = DenseVector[Double]((phikm1 - reverse(phikm1) * zeta(k)).toArray :+ zeta(k))
    }
    phik
  }

  def TacvfMA(theta: DenseVector[Double], maxlag: Int = 20): DenseVector[Double] = {
    if (theta.length == 0) {
      if (maxlag >= 0)
        return (DenseVector(1.0 +: Array.fill[Double](maxlag + 1)(0.0)))
      else sys.error("maxlag invalid")
    }
    val maxlagp1 = maxlag + 1
    val g = DenseVector.zeros[Double](maxlagp1)
    val th = DenseVector.vertcat(DenseVector(-1.0), theta)
    val qth = th.length
    val x = DenseVector.vertcat(th, DenseVector.zeros[Double](qth))
    val A = DenseMatrix.zeros[Double](qth, qth)
    val B = new DenseMatrix[Double](qth, qth, (abs(col(A) + row(A))).map(v => x(v)).toArray)
    val g1 = (B * th).toDenseVector
    if (g1.length < maxlagp1) g(0 until qth) := g1
    else g(0 until maxlagp1) := g1(0 until maxlagp1)
    g
  }

  def TacvfARMA(phi: DenseVector[Double], theta: DenseVector[Double], maxlag: Int = 20): DenseVector[Double] = {
    if (!(InvertibleQ(phi) & InvertibleQ(theta))) {
      System.err.println("TacvfARMA: Model is non-causal or non-invertible")
      System.err.println("phi =" + phi + ", theta = " + theta)
      return (null)
    }
    val p = phi.length
    val q = theta.length
    val maxlagp1 = maxlag + 1
    var res = DenseVector[Double](1.0 +: Array.fill(maxlagp1)(0.0))
    if (max(p, q) == 0) {
      return (res)
    }
    res = DenseVector.zeros[Double](maxlagp1)
    val r = max(p, q) + 1
    val b = DenseVector.zeros[Double](r)
    val C = DenseVector.zeros[Double](q + 1)
    C(0) = 1
    val theta2 = DenseVector[Double](-1.0 +: theta.toArray)
    val phi2 = DenseVector.zeros[Double](3 * r)
    phi2(r - 1) = -1
    if (p > 0) {
      phi2(r to r + p - 1) := phi
    }
    if (q > 0) {
      for (k <- 0 to q - 1) {
        C(k + 1) = -theta(k)
        if (p > 0) {
          for (i <- 0 to min(p, k + 1) - 1) {
            C(k + 1) = C(k + 1) + phi(i) * C(k - i)
          }
        }
      }
    }

    for (k <- 0 to q) {
      for (i <- k to q) {
        b(k) = b(k) - theta2(i) * C(i - k)
      }
    }

    if (p == 0) {
      res(0 to b.length - 1) := b
    } else {
      val a = DenseMatrix.zeros[Double](r, r)
      for (i <- 1 to r) {
        for (j <- 1 to r) {
          if (j == 1) {
            a(i - 1, j - 1) = phi2(r + i - 2)
          } else {
            a(i - 1, j - 1) = phi2(r + i - j - 1) + phi2(r + i + j - 3)
          }
        }
      }
      val g = a \ (-b)
      if (g.length <= maxlag) {
        res(0 to g.length - 1) := g
        for (i <- (r + 1) to maxlagp1) {
          res(i - 1) = phi dot reverse(res(i - p - 1 to i - 2))
        }
      } else if (g.length >= maxlagp1) {
        res = g(0 to maxlagp1 - 1)
      }
    }

    res
  }

  def ARToPacf(phi: DenseVector[Double]): DenseVector[Double] = {
    var phik = phi

    val L = phi.length
    if (L == 0)
      return (DenseVector.zeros[Double](0))
    val pi = DenseVector.zeros[Double](L)
    for (k <- 1 to L) {
      val LL = L + 1 - k
      val a = phik(LL - 1)
      pi(L - k) = a
      val phikp1 = DenseVector.vertcat(phik(0 until LL - 1), phik(LL to phik.length - 1))
      if (abs(a) == 1)
        System.err.println("transformation is not defined, partial correlation = 1")
      phik = (phikp1 + reverse(phikp1) * a) / (1 - pow(a, 2))
    }
    pi
  }

  def InvertibleQ(phi: DenseVector[Double]): Boolean = abs(ARToPacf(phi)).forall(v => v < 1)

  def PacfDL(c: DenseVector[Double], LinearPredictor: Boolean = false): Either[DenseVector[Double], (DenseVector[Double], DenseVector[Double], Double)] = {
    val L = c.length - 1
    val d = c(1 to L)
    var phik = DenseVector.zeros[Double](0)
    var vk = c(0)
    var pi = DenseVector.zeros[Double](0)
    if (L != 0) {
      phik = DenseVector[Double](Array(c(1) / c(0)))
      pi = DenseVector.zeros[Double](L)
      pi(0) = phik(0)
      vk = c(0) * (1 - pow(pi(0), 2))
    }

    if (L > 1) {
      for (k <- 2 to L) {
        val vkm1 = vk
        val phikm1 = phik
        val a = (DenseVector[Double](1.0 +: (-phikm1).toArray) dot reverse(d(0 to k - 1))) / vk
        phik = DenseVector[Double]((phikm1 - reverse(phikm1) * a).toArray :+ a)
        vk = vkm1 * (1 - pow(a, 2))
        pi(k - 1) = a
      }
    }
    if (!LinearPredictor)
      Left(pi)
    else Right((pi, phik, vk))
  }

  def FastLoglikelihoodAR(phi: DenseVector[Double], n: Int, CD: DenseMatrix[Double]): Double = {
    val phis = DenseVector[Double](1.0 +: (-phi).toArray)
    var LL = 10E-35
    try {
      LL = -log(DetAR(phi)) / 2.0 - (n / 2.0) * log(((phis.t * CD).t dot phis) / n)
    } catch {
      case e: Throwable => sys.error("error computing LL")
    }
    LL
  }

  def DetAR(phi: DenseVector[Double]): Double = {
    val z = ARToPacf(phi)
    val temp = -pow(z, 2) + 1.0
    1 / 1.to(temp.length).map(i => pow(temp(i - 1), i)).fold(1.0)((a, b) => a * b)
  }

  def gradientApproximate(f: DenseVector[Double] => Double, x: DenseVector[Double], h: Double): DenseVector[Double] = {
    val n = x.length
    val g = DenseVector.zeros[Double](n)
    for (i <- 0 until n) {
      val y = x.copy
      val z = x.copy
      y(i) += h
      g(i) = (f(y) - f(z)) / h
    }
    g
  }

  def gradientCenterApproximate(f: DenseVector[Double] => Double, x: DenseVector[Double], h: Double): DenseVector[Double] = {
    val n = x.length
    val g = DenseVector.zeros[Double](n)
    for (i <- 0 until n) {
      val y = x.copy
      val z = x.copy
      y(i) += h
      z(i) -= h
      g(i) = (f(y) - f(z)) / (2 * h)
    }
    g
  }

  def gradientCenter4Approximate(f: DenseVector[Double] => Double, x: DenseVector[Double], h: Double): DenseVector[Double] = {
    val n = x.length
    val g = DenseVector.zeros[Double](n)
    for (i <- 0 until n) {
      val y1 = x.copy
      val y2 = x.copy
      val y3 = x.copy
      val y4 = x.copy
      y1(i) += 2 * h
      y2(i) += h
      y3(i) -= h
      y4(i) -= 2 * h
      g(i) = (-f(y1) + 8 * f(y2) - 8 * f(y3) + f(y4)) / (12 * h)
    }
    g
  }

  def diff(x: DenseVector[Double], lag: Int = 1, differences: Int = 1): DenseVector[Double] = {
    val xlen = x.length
    if (lag < 1 || differences < 1) sys.error("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen) return (DenseVector.zeros[Double](0))
    var y = x.copy
    for (i <- 1 to differences) {
      val i1 = 0 until y.length - lag
      val i2 = lag until y.length
      y = y(i2) - y(i1)
    }
    y
  }

  def trunc(x: Double): Double = x.toInt.toDouble

  def floor(x: DenseVector[Double]): DenseVector[Double] = x.map(v => trunc(v))

  def row(x: DenseMatrix[Double]): DenseMatrix[Int] = {
    val rowMatrix = DenseMatrix.zeros[Int](rows = x.rows, cols = x.cols)
    for (i <- 0 until x.rows) {
      for (j <- 0 until x.cols) {
        rowMatrix(i, j) = i
      }
    }
    rowMatrix
  }

  def col(x: DenseMatrix[Double]): DenseMatrix[Int] = {
    val rowMatrix = DenseMatrix.zeros[Int](rows = x.rows, cols = x.cols)
    for (i <- 0 until x.rows) {
      for (j <- 0 until x.cols) {
        rowMatrix(i, j) = j
      }
    }
    rowMatrix
  }

  def lower_tri(x: DenseMatrix[Double], diag: Boolean = true): DenseMatrix[Boolean] = {
    if (diag)
      (row(x) - col(x)).map(v => v >= 0)
    else (row(x) - col(x)).map(v => v > 0)
  }

  def FromSymmetricStorageUpper(x: DenseVector[Double]): DenseMatrix[Double] = {
    val l = x.length
    val n = ((-1 + sqrt(1 + 8 * l)) / 2).toInt
    val z = lower_tri(DenseMatrix.zeros[Double](n, n), diag = true).map(v => if (v) 1.0 else 0.0)
    var count = 0
    for (i <- 0 until z.cols) {
      for (j <- 0 until z.rows) {
        if (z(j, i) != 0.0) {
          z(j, i) = x(count)
          count += 1
        }
      }
    }
    if (count != l - 1) println("warning: number of items to replace is not a multiple of replacement length")
    val ztranspose = z.t
    z + ztranspose - diag(diag(ztranspose))
  }

  def GetB(phi: DenseVector[Double]): DenseMatrix[Double] = {
    val p = phi.length
    var a = pow(phi, 2)
    val phistat = new StatDenseVector(phi)
    if (p != 1) {
      a = phistat.acf(Type = "covariance", demean = false) * p.toDouble
      for (i <- 1 until p) {
        val phistat2 = new StatDenseVector(phi(i until p).toArray ++ Array.fill[Double](i)(0))
        a = DenseVector.vertcat(a, (phistat2.acf(Type = "covariance", demean = false) * p.toDouble).apply(0 to p - i - 1))
      }
    }
    FromSymmetricStorageUpper(x = a)
  }

  def cumsum(x: DenseVector[Double]): DenseVector[Double] = {
    val cumSum = DenseVector.zeros[Double](x.length)
    var total = 0.0
    for (i <- 0 until x.length) {
      total += x(i)
      cumSum(i) = total
    }
    cumSum
  }

  def GetKappa(phi: DenseVector[Double]): DenseVector[Double] = {
    val tacvfMA = TacvfMA(theta = phi, maxlag = phi.length)
    reverse(cumsum(reverse(tacvfMA(1 until tacvfMA.length))))
  }

  def rowSums(x: DenseMatrix[Double]): DenseVector[Double] = {
    val rowSum = DenseVector.zeros[Double](x.cols)
    for (i <- 0 until x.cols) rowSum(i) = sum(x(::, i))
    rowSum
  }

  def Get1G(phi: DenseVector[Double], n: Int): DenseVector[Double] = {
    val p = phi.length
    val x0 = pow(sum(DenseVector[Double]((1.0 +: (-phi).toArray))), 2)
    val x = -rowSums(x = GetB(phi = phi)) - GetKappa(phi = phi) + x0
    DenseVector.vertcat(x, DenseVector.fill(n - 2 * p) { x0 }, reverse(x))
  }

  def GetARMeanMLE(z: DenseVector[Double], phi: DenseVector[Double]): Double = {
    if (z.length < 2 * phi.length) sys.error("the length of phi is coefficient of AR is too large, 2 * phi.length > data.length")
    val g1 = Get1G(phi = phi, n = z.length)
    (g1 dot z) / sum(g1)
  }

  def BackcastResidualsAR(y: DenseVector[Double], phi: DenseVector[Double], Q: Int = 100, demean: Boolean = true): DenseVector[Double] = {
    var z = y
    if (demean) z = y - mean(y)
    val p = phi.length
    var a = z
    if (p != 0) {
      val n = z.length
      val nQ = n + Q
      val zR = DenseVector.zeros[Double](nQ)
      val zF = DenseVector.zeros[Double](nQ)
      val e = DenseVector.zeros[Double](nQ)
      a = DenseVector.zeros[Double](nQ)
      val r = p + 1
      zR(0 until n) := reverse(z)
      zF(Q until Q + n) := z
      for (i <- r to n) e(i - 1) = zR(i - 1) - (phi dot reverse(zR(i - p - 1 to i - 2)))
      for (i <- 1 to Q) zR(n + i - 1) = phi dot reverse(zR(n + i - p - 1 to n + i - 2))
      zF(0 until Q) := reverse(zR(n to n + Q - 1))
      for (i <- r to nQ) a(i - 1) = zF(i - 1) - (phi dot reverse(zF(i - p - 1 to i - 2)))
      zF(Q until Q + n) := z
      a = a(Q until Q + n)
    }
    a
  }

  def LjungBoxTest(res: DenseVector[Double], k: Int = 0, maxlag: Int = 30, StartLag: Int = 1, SquaredQ: Boolean = false): DenseMatrix[Double] = {
    if (!(k >= 0 && StartLag >= 1 && maxlag >= StartLag)) sys.error("Should be k >= 0 and StartLag >= 1 and maxlag >= StartLag, please check again")
    val n = res.length
    val L0 = StartLag
    var z = res
    var kpar = k
    if (SquaredQ) {
      z = pow((res - mean(res)), 2)
      kpar = 0
    }
    val zacf = new StatDenseVector(z)
    val ra = zacf.acf(maxlag = maxlag)(1 to maxlag)
    val lags = DenseVector.rangeD(L0, maxlag + 1)
    val QQ = cumsum(pow(ra, 2) :/ reverse(DenseVector.rangeD(n - maxlag, n)))(L0 - 1 to maxlag - 1) * n.toDouble * (n + 2.0)
    val df = (lags - kpar.toDouble) map (v => if (v > 0) v else 1)
    val dn = df.length
    val pv = DenseVector.zeros[Double](dn)
    for (i <- 0 to dn - 1) {
      pv(i) = 1 - new ChiSquaredDistribution(df(i)).cumulativeProbability(QQ(i))
    }
    val a = DenseMatrix.zeros[Double](rows = dn, cols = 3)
    a(::, 0) := lags
    a(::, 1) := QQ
    a(::, 2) := pv
    a
  }

  def BoxCox(x: DenseVector[Double], lambda: Double): DenseVector[Double] = {
    var out = x.copy
    if (lambda < 0) out = out.map(v => if (v < 0) NaN else v)
    if (lambda == 0) out = log(out)
    else out = ((signum(out) :* pow(abs(out), lambda)) - 1.0) / lambda
    out
  }

  def InvBoxCox(x: DenseVector[Double], lambda: Double): DenseVector[Double] = {
    var out = x.copy
    if (lambda < 0) out = out.map(v => if (v > -1 / lambda) NaN else v)
    if (lambda == 0) out = exp(out)
    else {
      val xx = out * lambda + 1.0
      out = signum(xx) :* pow(abs(xx), (1 / lambda))
    }
    out
  }

  def embed(x: DenseVector[Double], dimension: Int = 1) = {
    val n = x.length
    if ((dimension < 1) | (dimension > n)) sys.error("wrong embedding dimension")
    val m = n - dimension + 1
    val data = DenseMatrix.zeros[Double](m, dimension)
    for (i <- 0 to m - 1) for (j <- 0 to dimension - 1) data(i, j) = x(dimension + i - j - 1)
    data
  }

  def adfTest(x: DenseVector[Double], alternative: String = "Stationary", lag: Int = 0): (Double, Double, String, Double) = {
    var k = lag
    if (lag == 0) k = pow((x.length - 1), (1.0 / 3)).toInt
    k += 1
    val y = diff(x)
    val n = y.length
    val z = embed(y, k)
    val yt = z(::, 0)
    val xt1 = x(k - 1 until n)
    val tt = DenseVector.rangeD(k, n + 1)
    val lm = new OLSMultipleLinearRegression()
    var xtt = xt1.toArray.zip(tt.toArray).map(v => Array(v._1, v._2))
    if (k > 1) {
      val yt1 = z(::, 1 until k).copy
      for (i <- 0 until k - 1) xtt = xtt.zip(yt1(::, i).toArray).map(v => v._1 :+ v._2)
    }
    lm.newSampleData(yt.toArray, xtt)
    val STAT = lm.estimateRegressionParameters()(1) / lm.estimateRegressionParametersStandardErrors()(1)
    val tablepositive = DenseMatrix((4.38, 4.15, 4.04, 3.99, 3.98, 3.96), (3.95,
      3.8, 3.73, 3.69, 3.68, 3.66), (3.6, 3.5, 3.45, 3.43,
      3.42, 3.41), (3.24, 3.18, 3.15, 3.13, 3.13, 3.12), (1.14,
      1.19, 1.22, 1.23, 1.24, 1.25), (0.8, 0.87, 0.9, 0.92,
      0.93, 0.94), (0.5, 0.58, 0.62, 0.64, 0.65, 0.66), (0.15,
      0.24, 0.28, 0.31, 0.32, 0.33))
    val table = -tablepositive.t
    val tablen = table.cols
    val tableT = DenseVector[Double](25, 50, 100, 250, 500, 1e+05)
    val test = DenseVector(4.38, 4.15, 4.04, 3.99, 3.98, 3.96)
    val tablep = DenseVector[Double](0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
    val tableipl = DenseVector.zeros[Double](tablen)
    for (i <- (0 until tablen)) {
      tableipl(i) = LinearInterpolator(tableT, table(::, i)).apply(n)
    }
    val interpol = LinearInterpolator(tableipl, tablep).apply(STAT)
    println("warning : p-value greater than printed p-value")
    var PVAL = interpol
    if (alternative == "explosive") PVAL = 1 - interpol
    val PARAMETER = k - 1
    (STAT, PARAMETER, alternative, PVAL)
  }

  def ma_inf(phi: DenseVector[Double], theta: DenseVector[Double], n: Int = 50): DenseVector[Double] = {
    if (n == 0) return DenseVector.ones[Double](1)
    val theta2 = DenseVector.vertcat(theta, DenseVector.zeros[Double](n))
    val p = phi.length
    val psi = DenseVector[Double]((Array.fill[Double](p)(0) :+ 1.0) ++ Array.fill[Double](n)(0))
    for (j <- 1 to n) psi(j + p) = theta2(j - 1) + (phi dot reverse(psi(j to p + j - 1)))
    psi(p to p + n)
  }

  def aacvf(phi: DenseVector[Double], theta: DenseVector[Double], sigma2: Double, h: Int): DenseVector[Double] = {
    val p = phi.length
    val q = theta.length
    val psi = ma_inf(phi, theta, q)
    val theta2 = DenseVector[Double](1.0 +: theta.toArray)
    def f1(k: Int): Double = theta2(k - 1 to q) dot psi(0 to (q - k + 1))
    val r = DenseVector.zeros[Double](max(p + 1, q + 1, h + 1))
    for (i <- 1 to q + 1) r(i - 1) = sigma2 * f1(i)
    val Ap = DenseMatrix.zeros[Double](p + 1, 2 * p + 1)
    def f2(k: Int): DenseVector[Double] = DenseVector[Double]((Array.fill[Double](p - k)(0) :+ 1.0) ++ (-phi).toArray ++ Array.fill[Double](k)(0))
    for (i <- 0 to p) Ap(i, ::) := f2(i).t
    val A = DenseMatrix.zeros[Double](p + 1, p + 1)
    A(::, 0) := Ap(::, p)
    A(::, 1 to p) := Ap(::, p + 1 to 2 * p) + Ap(::, p - 1 to 0 by -1)
    var gamma = DenseVector.zeros[Double](max(p + 1, h + 1))
    gamma(0 to p) := A \ r(0 to p)
    if (h > p) for (k <- p + 1 to h) gamma(k) = r(k) + (phi dot reverse(gamma(k - p to k - 1)))
    else if (h < p) gamma = gamma(0 to h)
    gamma
  }

  def innovations(x: DenseVector[Double], model: (DenseVector[Double], DenseVector[Double], Double)): DenseVector[Double] = {
    val (xhat, v) = innovation_kernel(x = x, model = model)
    x - xhat
  }

  def innovation_update(x: DenseVector[Double], model: (DenseVector[Double], DenseVector[Double], Double)): (Double, Double) = {
    val (xhat, v) = innovation_kernel(x = x, model = model)
    val n = x.length.toDouble 
    val sigma2 = sum(pow((x - xhat), 2) :/ v) / n
    val phi = model._1
    val theta = model._2
    var p = phi.length
    var q = theta.length
    val loglike = -(n / 2) * log(2 * constants.Pi * sigma2) - sum(log(v)) / 2 - n / 2
    val aicc = -2 * loglike + 2 * (p + q + 1) * n / (n - p - q - 2)
    (sigma2, aicc)
  }

  def innovation_kernel(x: DenseVector[Double], model: (DenseVector[Double], DenseVector[Double], Double)) = {
    val phi = model._1
    val theta = model._2
    val sigma2 = model._3
    val N = x.length
    val theta_r = DenseVector[Double]((1.0 +: theta.toArray) ++ Array.fill[Double](N)(0))
    val gamma = aacvf(phi, theta, sigma2, N - 1)
    val p = phi.length
    val q = theta.length
    val m = max(p, q)

    def kappa(i: Int, j: Int): Double = {
      if (j > m) return theta_r(0 to q) dot theta_r(i - j to i - j + q)
      else if (i > 2 * m) return 0
      else if (i > m) {
        var sum = 0.0
        if (p >= 1) {
          for (k <- 1 - i + j to p - i + j) {

            sum += phi(k - 1 + i - j) * gamma(abs(k))
          }
        } else {
          for (k <- p - i + j to 1 - i + j) {
            sum += phi(k - 1 + i - j) * gamma(abs(k))
          }
        }
        return (gamma(i - j) - sum) / sigma2
      } else return gamma(i - j) / sigma2
    }
    val Theta = DenseMatrix.zeros[Double](N - 1, N - 1)
    val v = DenseVector.zeros[Double](N)
    v(0) = kappa(1, 1)
    for (n <- 1 until N) {
      for (k <- 0 until n) {
        val u = kappa(n + 1, k + 1)
        var sum = 0.0
        if (k > 0) {
          for (j <- 0 to k - 1) {
            sum += Theta(k - 1, k - j - 1) * Theta(n - 1, n - j - 1) * v(j)
          }
        }

        Theta(n - 1, n - k - 1) = (u - sum) / v(k)
      }
      val s = reverse(pow(Theta(n - 1, 0 to n - 1).t, 2)) dot v(0 to n - 1)
      v(n) = kappa(n + 1, n + 1) - s
    }

    val xhat = DenseVector.zeros[Double](N)
    if (m > 1) {
      for (n <- 1 to m - 1) {
        xhat(n) = Theta(n - 1, 0 until n).t dot reverse(x(0 until n) - xhat(0 until n))
      }
    }

    for (n <- m to N - 1) {
      val A = phi dot reverse(x(n - p to n - 1))
      var B = 0.0
      for (i <- 1 to q) {
        B += Theta(n - 1, i - 1) * (x(n - i) - xhat(n - i))
      }
      xhat(n) = A + B
    }
    (xhat, v)
  }

  def forecast_arma(y: DenseVector[Double], model: (DenseVector[Double], DenseVector[Double], Double), h: Int, demean : Boolean = true): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = y.length
    var mu = mean(y)
    if (!demean) mu = 0.0
    var x = y - mu
    val phi = model._1
    val theta = model._2
    var dx = innovations(x = x, model = model)
    dx = DenseVector.vertcat(dx, DenseVector.zeros[Double](h))
    x = DenseVector.vertcat(x, DenseVector.zeros[Double](h))
    val p = phi.length
    val q = theta.length
    for (t <- n to n + h - 1) {
      val A = phi dot reverse(x(t - p to t - 1))
      val B = theta dot reverse(dx(t - q to t - 1))
      x(t) = A + B
    }
    val pred = x(n to n + h - 1) + mu
    (pred, phi, null, null)
  }
  
  def forecast_diff(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int)
   : (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = x.length
    val lag = xv(k).toInt
    val differences = xv(k + 1).toInt
    val y = diff(x, lag, differences)
    var (pred, phi, l, u) = forecast_transform(x = y, xv = xv, model = model, h = h, k = k + 3)
    pred = diffinv(pred, lag, differences, x(n - lag to n - 1))
    pred = pred(lag to lag + h - 1)
    if (phi == null) sys.error("diff before log")
    phi = DenseVector[Double](1.0 +: (-phi).toArray)
    phi = DenseVector.vertcat(phi, DenseVector.fill[Double](lag, 0.0)) - DenseVector.vertcat(DenseVector.fill[Double](lag, 0.0), phi)
    phi = -phi(1 until phi.length)
    (pred, phi, l, u)
  }

  def forecast_log(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int)
   : (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    var (pred, phi, l, u) = forecast_transform(x = log(x), xv = xv, model = model, h = h, k = k + 1)
    if (phi != null) {
      val theta = model._2
      val sigma2 = model._3
      val psi = ma_inf(phi = phi, theta = theta, n = h)
      def g(j: Int) = sum(pow(psi(0 to j - 1), 2))
      val se = DenseVector.zeros[Double](h)
      for (i <- 1 to h) se(i - 1) = sqrt(sigma2 * g(i))
      l = pred - (se * 1.96)
      u = pred + (se * 1.96)
    }
    pred = exp(pred)
    l = exp(l)
    u = exp(u)
    (pred, null, l, u)
  }
  
  def forecast_season(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int = 1)
  	: (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = x.length	
    val d = xv(k).toInt
    val q = d / 2
    def F1(t: Int): Double = x(t - q - 1) / 2 + sum(x(t - q to t + q - 2)) + x(t + q - 1) / 2
    var m = DenseVector.zeros[Double](n - 2 * q)
    for (i <- q + 1 to n - q) m(i - q - 1) = F1(i) / d
    m = DenseVector[Double](Array.fill[Double](q)(0) ++ m.toArray ++ Array.fill[Double](q)(0))
    if (d != 2 * q) m = smooth_ma(x, q)
    val dx = x - m
    def F2(k: Int): Double = mean(dx(k + q - 1 to n - q - 1 by d))
    var w = DenseVector.zeros[Double](d)
    for (i <- 1 to d) w(i - 1) = F2(i)
    w = w - mean(w)
    var s = w
    for (i <- 1 until n + q + h) s = DenseVector.vertcat(s, w)
    s = s(q to q + n + h - 1)
    val y = x - s(0 until n)
    var (pred, phi, l, u) = forecast_transform(x = y, xv = xv, model = model, h = h, k = k + 2)
    pred = pred + s(n to n + h - 1)
    if (phi == null) {
        l = l + s(n to n + h - 1)
        u = u + s(n to n + h - 1)
    }
    (pred, phi, l, u)
  }
  
  def forecast_trend(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int = 1)
  	: (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = x.length
    val p = xv(k).toInt
    val X = DenseMatrix.zeros[Double](n + h, p + 1)
    for (j <- 0 to p) X(::, j) := pow(DenseVector.rangeD(1.0, n + h + 1), j)
    val b = X(0 until n, ::) \ x
    val xhat = X * b
    val y = x - xhat(0 until n)
    var (pred, phi, l, u) = forecast_transform(x = y, xv = xv, model = model, h = h, k = k + 2) 
    pred = pred + xhat(n until n + h)
    if (phi == null){
      l = l + xhat(n until n + h)
      u = u + xhat(n until n + h)
    }
    (pred, phi, l, u)
  }
  
  def forecast_transform(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int = 1)
  : (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    if (k > xv.length) return forecast_arma(y = x, model = model, h = h)
    if (xv(k - 1) == "diff") return forecast_diff(x = x, xv = xv, model = model, h = h, k = k)
    if (xv(k - 1) == "log") return forecast_log(x = x, xv = xv, model = model, h = h, k = k)
    if (xv(k - 1) == "season") return forecast_season(x = x, xv = xv, model = model, h = h, k = k)
    if (xv(k - 1) == "trend") return forecast_trend(x = x, xv = xv, model = model, h = h, k = k)
    else sys.error("x vector transformation is invalid")
  }

  def forecast(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int = 1)
  : (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
  	var (pred, phi, l, u) = forecast_transform(x = x, xv = xv, model = model, h = h, k = 1) 
  	if (phi != null){
  	  val theta = model._2
  	  val psi = ma_inf(phi = phi, theta = theta, n = h)
  	  val sigma2 = model._3
  	  def g(j: Int) : Double = sum(pow(psi(0 until j), 2))
  	  val se = DenseVector.zeros[Double](h)
  	  for (i <- 1 to h) se(i) = sqrt(sigma2 * g(i))
  	  l = pred - (se * 1.96)
  	  u = pred + (se * 1.96)
  	}
  	(pred, phi, l, u)
  }
    
  def season(x: DenseVector[Double], d: Int) : DenseVector[Double] = {
    val n = x.length
    val q = d / 2
    def F1(t: Int): Double = x(t - q - 1) / 2 + sum(x(t - q to t + q - 2)) + x(t + q - 1) / 2
    var m = DenseVector.zeros[Double](n - 2 * q)
    for (i <- q + 1 to n - q) m(i - q - 1) = F1(i) / d
    m = DenseVector[Double](Array.fill[Double](q)(0) ++ m.toArray ++ Array.fill[Double](q)(0))
    if (d != 2 * q) m = smooth_ma(x, q)
    val dx = x - m
    def F2(k: Int): Double = mean(dx(k + q - 1 to n - q - 1 by d))
    var w = DenseVector.zeros[Double](d)
    for (i <- 1 to d) w(i - 1) = F2(i)
    w = w - mean(w)
    var s = w
    for (i <- 1 until n + q) s = DenseVector.vertcat(s, w)
    s(q to q + n - 1)
  }

  def smooth_ma(x: DenseVector[Double], q: Int): DenseVector[Double] = {
    val n = x.length
    val y = DenseVector[Double](Array.fill[Double](q)(x(0)) ++ x.toArray ++ Array.fill[Double](q)(x(n - 1)))
    def F(t: Int): Double = sum(y(t - q - 1 to t + q - 1)) / (2 * q + 1)
    var m = DenseVector.zeros[Double](n)
    for (i <- q + 1 to n + q) m(i - q - 1) = F(i)
    m
  }

  def Resid(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double) = null): DenseVector[Double] = {
    var y = x
    var k = 1
    while (k < xv.length) {
      if (xv(k - 1) == "diff") {
        val lag = xv(k).toInt
        val differences = xv(k + 1).toInt
        y = diff(y, lag, differences)
        k += 3
      } else if (xv(k - 1) == "log") {
        y = log(y)
        k += 1
      } else if (xv(k - 1) == "season") {
        val d = xv(k).toInt
        y = y - season(y, d)
        k = k + 2
      } else if (xv(k - 1) == "trend") {
        val p = xv(k).toInt
        y = y - trend(y, p)
        k = k + 2
      } else {
        sys.error("x vector transformation is invalid")
      }
    }
    y = y - mean(y)
    if (model != null) {
      val (xhat, v) = innovation_kernel(y, model)
      y = (y - xhat) / sqrt(v)
    }
    y
  }

  def trend(x: DenseVector[Double], p: Int): DenseVector[Double] = {
    val n = x.length
    val X = DenseMatrix.zeros[Double](n, p + 1)
    for (i <- 0 to p) X(::, i) := pow(DenseVector.rangeD(1, n + 1), i)
    val b = X \ x
    val xhat = X * b
    xhat
  }

  def diffinv(x: DenseVector[Double], lag: Int = 1, differences: Int = 1, xi: DenseVector[Double] = null): DenseVector[Double] = {
    if (lag < 1 || differences < 1) sys.error("bad value for 'lag' or 'differences'")
    val n = x.length
    val y = DenseVector.zeros[Double](n + lag)
    if (xi.length != lag * differences) sys.error("'xi' has not the right length")
    if (differences == 1) {
      if (xi != null) {
        y(0 until lag) := xi(0 until lag)
      }
      for (i <- lag until n + lag) {
        y(i) = y(i - lag) + x(i - lag)
      }
      y
    } else {
      return diffinv(diffinv(x, lag, differences - 1, diff(xi, lag = lag, differences = 1)), lag, 1, xi(0 until lag))
    }
  }

  def predict(predictionLength: Int = 12): Array[Double] = {
    null
  }
  
  def main(args: Array[String]) {
    val x = csvReader("D://Data/sharefolder/testdata.csv")
    val y = DenseVector.rangeD(-10, 21)
    val zeta = DenseVector[Double](Array(0.4, 0.3, 0.2, 0.2, 0.3, 0.4))
    val phi = zeta(0 to 2)
    val theta = zeta(3 to 5)
    val features = DenseMatrix.create(9, 2, Array(1.0, 2, 5, 7, 9, 13, 15, 19, 23, 1, 2, 3, 4, 5, 6, 7, 8, 9))
    val test = DenseVector(1.0, 2, 5, 7, 9, 13, 15, 19, 23)
    val model = (phi, theta, 0.5)
    println(innovation_update(DenseVector(x),model))
//    println(forecast_transform(DenseVector(x), List("diff","1", "1"), model, 10, 1))
    /*    val CD = ChampernowneD(y, 20, true)
    val g = TacvfARMA(PacfToAR(phi), PacfToAR(theta),20)
    val xpar = PacfDL(g, true).right.get._2
    val xx = DenseVector(1.0,1.0)
    def f(x : DenseVector[Double]) : Double = sum(x)
    println(gradientApproximate(f,xx,0.1))
    
    val xs = DenseMatrix.zeros[Double](2,3)*/
    //    val res = BackcastResidualsAR(y, zeta)
    //   println(LjungBoxTest(res, maxlag = 10, k = 5))
  }
}