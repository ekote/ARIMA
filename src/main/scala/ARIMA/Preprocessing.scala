
/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

package ARIMA

import breeze.linalg.DenseVector
import stat.StatDenseVector
import breeze.numerics.sqrt
import breeze.stats.mean
/*import breeze.stats._
import org.apache.commons.math3.distribution.ChiSquaredDistribution
*/


abstract class Preprocessing(data: Array[Double]) {

  def this(y: StatDenseVector) {
    this(y.toArray)
  }

  val x = new DenseVector(data)

  // time series mean
  val mu: Double = mean(x)

  // mean zero time series
  val y = new StatDenseVector(x - mu toArray)

  // the size of Time series data  
  val n: Int = x.length

  // Maximum lag to consider  
  val maxOrder: Int = scala.math.min((20. * math.log10(n)).toInt, n - 1)

  // Time series variance
  val sig2: Double = y.sig

  // Time series standard deviation
  val stddev: Double = sqrt(sig2)

  // Auto-covariance vector, c(k) = cov[X0,Xk]
  val acvf: DenseVector[Double] = y.acf(Type = "covariance")

  // Auto-Correlation Function (ACF)
  val acf: DenseVector[Double] = y.acf(Type = "correlation")

  // Partial Auto-Correlation Function (PACF)
  var pacf: DenseVector[Double] = DenseVector.zeros[Double](maxOrder + 1)
  
}