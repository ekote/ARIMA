package stat

import breeze.linalg.DenseVector

// you may use gradience approximation here to elevate the accuracy, default is using breeze package
object GradientApproximation {
    
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
}