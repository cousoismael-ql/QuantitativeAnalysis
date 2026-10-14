# Calibracion curva OIS (CLP Camara)

Construccion de una curva de descuento y proyeccion a partir de swaps OIS camara en pesos chilenos.

Incluye:
- Calendario de dias habiles con feriados de Chile y fechas de inicio T+2.
- Generacion de flujos (cupon cero hasta 18M, semestral desde 2Y), act/360.
- Calibracion por minimos cuadrados (`scipy.optimize.least_squares`) de zero rates para curva de descuento y proyeccion, con interpolacion spline sobre log(DF).
- Bootstrap secuencial de factores de descuento como comparacion.
- Forwards diarios e interpolacion log-lineal de DF.

## Datos

Los archivos `mkt_clpcam.csv` y `mkt_sofr.csv` contienen tasas **ficticias**, solo para que el notebook corra. No corresponden a cotizaciones reales de mercado.

`CLholidays.csv` se genero con la libreria [holidays](https://pypi.org/project/holidays/) (feriados publicos y bancarios de Chile, solo dias habiles).

## Uso

```
pip install -r requirements.txt
jupyter notebook ois_curve_calibrator.ipynb
```
