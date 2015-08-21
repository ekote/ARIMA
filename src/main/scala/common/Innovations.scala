/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */


package common

import breeze.linalg._
import breeze.numerics.{ pow, log, constants, abs }
import common.ARIMAUtils._

// this object is to provide Innovation algorithm
object Innovations {

  def innovations(x: DenseVector[Double], model: (DenseVector[Double], DenseVector[Double], Double)): DenseVector[Double] = {
    val (xhat, v) = innovation_kernel(x = x, model = model)
    x - xhat
  }

  // calculate aicc with innovations algorithm
  def innovation_update(x: DenseVector[Double], model: (DenseVector[Double], DenseVector[Double], Double)): (Double, Double) = {  
    val model2 = model.copy( _3 = 1.0)
    val (xhat, v) = innovation_kernel(x = x, model = model)
    val n = x.length.toDouble 
    val sigma2 = sum(pow((x - xhat), 2) :/ v) / n
    val phi = model._1
    val theta = model._2
    var p = phi.length
 //   if (!any(phi)) p = 0
    var q = theta.length
 //   if (!any(theta)) q = 0
    val loglike = -(n / 2) * log(2 * constants.Pi * sigma2) - sum(log(v)) / 2 - n / 2
    val aicc = -2 * loglike + 2 * (p + q + 1) * n / (n - p - q - 2)
    (sigma2, aicc)
  }
  
  // calculate kernel in innovations algorithm
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
        for (k <- 1 - i + j to p - i + j) sum += phi(k - 1 + i - j) * gamma(abs(k))
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
      var A = 0.0
      if (phi.length != 0) A = phi dot reverse(x(n - p to n - 1))
      var B = 0.0
      for (i <- 1 to q) {
        B += Theta(n - 1, i - 1) * (x(n - i) - xhat(n - i))
      }
      xhat(n) = A + B
    }
    (xhat, v)
  }

}