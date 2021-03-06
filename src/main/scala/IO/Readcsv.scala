/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */
package IO

import scala.io.{Source}

// a simple csv reader
object Readcsv {

  def csvReader(path: String) = {
    //	val src = Source.fromFile("D://workspace/ARIMA/src/test/resources/timeseries_ppi.csv")
    //  val list = src.getLines().map(_.split(",")).toList.drop(1)
    val src = Source.fromFile(path)
    val list = src.getLines().map(_.split(",")).toList.drop(1)
    val data = list.map(x => x(1).toDouble).toArray
    
    data
  }
}