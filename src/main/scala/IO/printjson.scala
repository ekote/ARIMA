/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */
package IO

import IO.ParseJson._
import java.io._
import IO.DateFrequency._

// this object is to parse Json and output into a file
object printjson {

  
  def printToFile(f: java.io.File)(op: java.io.PrintWriter => Unit) {
    val p = new java.io.PrintWriter(f)
    try { op(p) } finally { p.close() }
  }

  def main(args: Array[String]) {
    val data = parseJson(path = "D://sql/MADBCN.json", startdate = "2011-01-10", enddate = "2013-10-10")
    val x = data.map(v => v._2)
    val x_weekly = FromDaily(x, 7)
    val x_monthly = FromDaily(x, 30)
    printToFile(new File("D://sql/MADBCNtrainweekly2.txt")) { p =>
      p.println("numberofpax")
      x_weekly.foreach(p.println)
    }
  }
}