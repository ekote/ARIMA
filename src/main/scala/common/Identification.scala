package common

import breeze.linalg.{ DenseVector, DenseMatrix, sum, reverse }
import breeze.numerics.{ log, signum, pow, sqrt, abs, NaN, exp }
import breeze.stats.mean
import breeze.interpolation.LinearInterpolator
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression
import org.apache.commons.math3.distribution.ChiSquaredDistribution
import common.ARIMAUtils.{ embed, diff, smooth_ma, cumsum }
import common.Innovations._
import stat.StatDenseVector

object ARIMAIdTransUtils {

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
      val yt1 = z(::, 1 until z.cols)
      for (i <- 0 until yt1.cols) xtt = xtt.zip(yt1(::, i).toArray).map(v => v._1 :+ v._2)
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
    var interpol = LinearInterpolator(tableipl, tablep).apply(STAT)
    if (interpol <= 0) interpol = 0.01
    if (interpol > 1) interpol = 1.0
    println("warning : p-value greater than printed p-value")
    var PVAL = interpol
    if (alternative == "explosive") PVAL = 1 - interpol
    val PARAMETER = k - 1
    (STAT, PARAMETER, alternative, PVAL)
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

  def season(x: DenseVector[Double], d: Int): DenseVector[Double] = {
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

  def trend(x: DenseVector[Double], p: Int): DenseVector[Double] = {
    val n = x.length
    val X = DenseMatrix.zeros[Double](n, p + 1)
    for (i <- 0 to p) X(::, i) := pow(DenseVector.rangeD(1, n + 1), i)
    val b = X \ x
    val xhat = X * b
    xhat
  }

  def Resid(x: DenseVector[Double], d: Int = 0, demean: Boolean = true, xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double) = null): DenseVector[Double] = {
    var y = x
    var k = 1
    var differences = d
    while (k < xv.length) {
      if (differences != 0) {
        if (xv(k - 1) == "diff") {
          val lag = xv(k).toInt
          y = diff(y, lag = lag, differences = differences)
          differences = 0
          k += 2
        } else {
          y = diff(y, lag = 1, differences = differences)
          differences = 0
        }
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
        println(xv)
        sys.error("x vector transformation is invalid, please verify diff and d")
      }
    }
    if (demean) y = y - mean(y)
    if (model != null) {
      val (xhat, v) = innovation_kernel(y, model)
      y = (y - xhat) / sqrt(v)
    }
    y
  }
}