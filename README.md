# AdaptiveQuadrature
Adaptive Gauss-Kronrod Quadrature for Java:

```java
UnaryOperator<double[]> parabola = x -> {
	int n = x.length;
	double[] f = new double[n];
	for (int i=0; i<n; ++i) {
		f[i] = 3.0*x[i]*x[i];
	}
	return f;
};

double[] result = AdaptiveQuadrature.integrate(parabola, 0.0, 1.0, 1.0e-6, 1.0e-6, 0);
System.out.printf("\\int_0^1 3*x^2 dx = %.3e +/- %.3e\n", result[0], result[1]);
```

## Maven coordinates

```
<dependency>
	<groupId>de.labathome</groupId>
	<artifactId>AdaptiveQuadrature</artifactId>
	<version>1.0.1</version>
</dependency>
```
