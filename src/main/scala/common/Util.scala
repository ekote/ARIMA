package common

import breeze.linalg.{ DenseVector, DenseMatrix, reverse, max, min, sum, diag }
import breeze.stats.{ mean }
import breeze.numerics.{ abs, pow, sqrt, log }
import stat.StatDenseVector
import common.Innovations._

object ARIMAUtils {

  def trunc(x: Double): Double = x.toInt.toDouble

  def floor(x: DenseVector[Double]): DenseVector[Double] = x.map(v => trunc(v))

  def dot(v1: DenseVector[Double], v2: DenseVector[Double]): Double = {
    var sum = 0.0
    for (i <- 0 until v1.length) sum += v1(i) * v2(i)
    sum
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
  
  def ARErrEst(y: DenseVector[Double], t: Int, phi: DenseVector[Double]): Double = y(t) - (phi dot reverse(y(t - phi.length to t - 1)))

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
    for (j <- 0 to p) for (i <- 0 to j - 1) E(j, i) = E(i, j)
    A - E
  }

  def PacfToAR(zeta: DenseVector[Double]): DenseVector[Double] = {
    val L = zeta.length
    if (L == 0) return (DenseVector.zeros[Double](0))
    if (L == 1) return (zeta)
    var phik = DenseVector[Double](Array(zeta(0)))
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
      println("TacvfARMA: Model is non-causal or non-invertible")
      sys.error("phi =" + phi + ", theta = " + theta)
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

  def DetAR(phi: DenseVector[Double]): Double = {
    val z = ARToPacf(phi)
    val temp = -pow(z, 2) + 1.0
    1 / 1.to(temp.length).map(i => pow(temp(i - 1), i)).fold(1.0)((a, b) => a * b)
  }

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
      case e: Throwable => e.getCause
    }
    LL
  }

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

  def embed(x: DenseVector[Double], dimension: Int = 1) = {
    val n = x.length
    if ((dimension < 1) | (dimension > n)) sys.error("wrong embedding dimension")
    val m = n - dimension + 1
    val data = DenseMatrix.zeros[Double](m, dimension)
    for (i <- 0 to m - 1) for (j <- 0 to dimension - 1) data(i, j) = x(dimension + i - j - 1)
    data
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

  def ma_inf(phi: DenseVector[Double], theta: DenseVector[Double], n: Int = 50): DenseVector[Double] = {
    if (n == 0) return DenseVector.ones[Double](1)
    val theta2 = DenseVector.vertcat(theta, DenseVector.zeros[Double](n))
    val p = phi.length
    val psi = DenseVector[Double]((Array.fill[Double](p)(0) :+ 1.0) ++ Array.fill[Double](n)(0))
    for (j <- 1 to n) psi(j + p) = theta2(j - 1) + (phi dot reverse(psi(j to p + j - 1)))
    psi(p to p + n)
  }

  def smooth_ma(x: DenseVector[Double], q: Int): DenseVector[Double] = {
    val n = x.length
    val y = DenseVector[Double](Array.fill[Double](q)(x(0)) ++ x.toArray ++ Array.fill[Double](q)(x(n - 1)))
    def F(t: Int): Double = sum(y(t - q - 1 to t + q - 1)) / (2 * q + 1)
    var m = DenseVector.zeros[Double](n)
    for (i <- q + 1 to n + q) m(i - q - 1) = F(i)
    m
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



}