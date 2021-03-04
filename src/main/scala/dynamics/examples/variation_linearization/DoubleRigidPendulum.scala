package dynamics.examples.variation_linearization

import dynamics.calculus.MatrixManipulation.extractVariationCoefficients
import dynamics.coder.Matlab
import dynamics.coder.Latex.{print2LatexFile, printVLatex, variationCoeffs2LatexEquation, saveString2File}
import dynamics.data_types._

object DoubleRigidPendulum {
  def main() {

    //     RigidPendulum
    val filename: String = "VariationDoubleRigidPendulum"

    // define constant scalars
    val g = ConstScalar("g") // gravitational constant
    val m1 = ConstScalar("m1") // mass of pendulum 1
    val m2 = ConstScalar("m2") // mass of pendulum 2
    val mt = ConstScalar("mt")

    // define constant vectors
    val e3 = ConstVector("e3")
    val l1 = ConstVector("l1")
    val lc1 = ConstVector("lc1")
    val l2 = ConstVector("l2")
    val lc2 = ConstVector("lc2")

    // define constant matrices
    val J1 = ConstMatrix("J1")
    val J2 = ConstMatrix("J2")
    val I = IdentityMatrix()

    // define states
    val R1 = SO3("R1")
    val R2 = SO3("R2")
    val Om1 = R1.getTangentVector
    val Om2 = R2.getTangentVector

    val eta1 = R1.getVariationVector
    val eta2 = R2.getVariationVector
    val u = Vector("u")

    // set configuration variables
    val configVars = Tuple3(List(), List(), List(R1, R2))

    val d11 = J1 - CrossMap(lc1) *** CrossMap(lc1) * m1 - CrossMap(l1) *** CrossMap(l1) * (m1 + mt)
    val d12 = CrossMap(l1) *** R1.T *** R2 *** CrossMap(lc2) * m2 * NumScalar(-1.0)
    val d21 = CrossMap(lc2) *** R2.T *** R1 *** CrossMap(l1) * m2 * NumScalar(-1.0)
    val d22 = J2 - CrossMap(lc2) *** CrossMap(lc2) * m2
    val D = List(List(d11, d12), List(d21, d22))

    val c1 = (CrossMap(Om1) *** d11 ** Om1) + CrossMap(l1) *** R1.T *** R2 *** CrossMap(Om2) *** CrossMap(Om2) ** lc2 * m2
    val c2 = (CrossMap(Om2) *** d22 ** Om2) + CrossMap(lc2) *** R2.T *** R1 *** CrossMap(Om1) *** CrossMap(Om1) ** l1 * m2
    val C = List(c1, c2)

    val g1 = CrossMap(lc1) *** R1.T ** e3 * Mul(m1, g) * NumScalar(-1.0) - CrossMap(l1) *** R1.T ** e3 * Mul(g, mt + m2)
    val g2 = CrossMap(lc2) *** R2.T ** e3 * Mul(m2, g)
    val G = List(g1, g2)

    val b1 = R1.T *** R2 * NumScalar(-1.0)
    val b2 = I
    val B = List(b1, b2)

    // equation for dOmega1
    val eq1 = d11 ** Om1.diff() + d12 ** Om2.diff() + c1 + g1 - b1 ** u
    // equation for dOmega2
    val eq2 = d21 ** Om1.diff() + d22 ** Om2.diff() + c2 + g2 - b2 ** u

    // list of variables to be extracted from the equation
    val variables = List(Om1.diff().delta(), Om2.diff().delta(), Om1.delta(), Om2.delta(), eta1, eta2, u.delta())

    var startTime = System.nanoTime() // track computation time

    val var_eqn1 = eq1.delta().basicSimplify()
    val var_eqn2 = eq2.delta().basicSimplify()
    val coeff1 = extractVariationCoefficients(var_eqn1, (List(), variables, List()))
    val coeff2 = extractVariationCoefficients(var_eqn2, (List(), variables, List()))

    val output_variables: List[Any] = List(Om1.diff().delta(), Om2.diff().delta())
    val state_variables: List[Any] = List(eta1, eta2, Om1.delta(), Om2.delta())
    val input_variables: List[Any] = List(u.delta())
    val dynamics = Matlab.generateLinearDynamics(List(coeff1, coeff2), output_variables, state_variables, input_variables)

    var endTime = System.nanoTime()
    println("ComputationTime: " + (endTime - startTime) / 1000000 + " ms")


    // Output
    var eqn_latex: String = variationCoeffs2LatexEquation(coeff1) + "\n\n" + variationCoeffs2LatexEquation(coeff2)


    saveString2File(eqn_latex, filename)
    Matlab.generateMatlabFunction(dynamics, filename)
    println("Testing done")


  }

}