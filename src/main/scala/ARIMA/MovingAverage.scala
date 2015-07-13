
/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

// package TimeSeriesForecast.ARIMA

package ARIMA

import breeze.linalg.{ DenseMatrix, DenseVector }
import breeze.numerics.pow
import stat.StatDenseVector

// Implementation of the Innovations Algorithm to estimate the coefficient of MA model.

class MovingAverage(data: Array[Double], q: Int, recursion_level: Int) extends Preprocessing(data: Array[Double]) {

  private val theta = DenseMatrix.zeros[Double](maxOrder + 1, maxOrder + 1)
  private val v = DenseVector.zeros[Double](maxOrder + 1)
  private val setheta = DenseVector.zeros[Double](q)

  def checkOrder: Boolean = q > maxOrder || q < 1

  def checkRecursionLevel: Boolean = recursion_level < q

  def Innovations_algorithm(data: Array[Double] = data, q: Int = q, recursion_level: Int = recursion_level): DenseVector[Double] = {

    v(0) = acvf(0)
    for (t <- 1 to maxOrder) {
      for (k <- 0 until t) {
        var sum = 0.0
        for (j <- 0 to k - 1) {
          sum += theta(k, k - j) * theta(t, t - j) * v(j)
        }
        theta(t, t - k) = (acvf(t - k) - sum) / v(k)
      }
      var sum2 = 0.0
      for (i <- 0 to t - 1) {
        sum2 += pow(theta(t, t - i), 2) * v(i)
      }
      v(t) = acvf(0) - sum2
    }
    setheta(0) = math.sqrt(1. / n)
    if (q > 1) {
      for (j <- 1 until q) {
        var sum3 = 0.0
        for (i <- 1 to j) {
          sum3 += pow(theta(recursion_level, i), 2)
        }
        setheta(j) = math.sqrt((1.0 + sum3) / n)
      }
    }
    theta(recursion_level, 1 to q).t
  }

  def run(): MovingAverageModel = {
    if (checkOrder) sys.error("Order p is too large, larger than max order min(n-1, 20*log10(n)) or less than 1")

    if (checkRecursionLevel) sys.error("Recursion Level is too small, less than q")

    val opttheta = Innovations_algorithm()

    val sigma2 = v(recursion_level)
    //   val Innovation = new Innovations(y,opttheta) not finished yet
    //   val (sigma2,aicc) = Innovation.update

    //   println("opttheta =" + opttheta.toArray.toList)
    //   println("sigma^2 = " + sigma2)
    //   println ("se_theta = " + setheta)

    new MovingAverageModel(y, mu, (opttheta, sigma2, setheta))

  }
}

class MovingAverageModel(y: StatDenseVector, mu: Double, result: (DenseVector[Double], Double, DenseVector[Double])) {

  val theta = result._1
  val sigma2 = result._2
  val se_theta = result._3

  def predict(predictionLength: Int = 12): Array[Double] = {
    null
  }

  def print_Result(): Unit = {
    println(result)
  }
}

object MovingAverage {

  def train(data: Array[Double], q: Int = 1, recursion_level: Int = 1): MovingAverageModel = {
    new MovingAverage(data, q, recursion_level).run()
  }
}