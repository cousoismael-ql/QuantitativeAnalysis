# Calibracion curvas OIS (CLP Camara / SOFR)

Construccion de curvas de descuento y proyeccion a partir de swaps OIS camara en pesos chilenos y SOFR en dolares, mas una curva basis entre ambas.

## Notebooks

**`double_currency_calibrator.ipynb`** - dos monedas (CLP Camara / SOFR)
- Curvas CLP Camara (semestral) y SOFR (anual) con sus calendarios.
- Swaps basis CLP vs SOFR y calibracion de la curva basis usando los forwards de ambas curvas.
- Graficos de factores de descuento y forwards SOFR.

## Datos

Todos los archivos de tasas son **ficticios**, solo para que los notebooks corran. No corresponden a cotizaciones reales de mercado.

- `mkt_clpcam.csv`, `mkt_sofr.csv`: tasas par por tenor.
- `curvasBasis.csv`: tenor, basis (pb), CLP Camara (%), SOFR (%), sin encabezado.
- `CLholidays.csv`, `USholidays.csv`: generados con la libreria [holidays](https://pypi.org/project/holidays/) (Chile publicos y bancarios, US federales), solo dias habiles.

## Uso

```
pip install -r requirements.txt
jupyter notebook
```
