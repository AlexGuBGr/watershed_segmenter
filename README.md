# watershed_segmenter
Small Python (v3.10.6) scripts to estimate the number and sizes of adipocytes in sets of adipocyte images. Python packages used for the image preparation and segmentation were OpenCV, Scikit-image, SciPy and Numpy. Python packages used for visualization and data compilation were Matplotlib, Seaborn and Pandas.  

## Run segmentation.py
```python segmentation.py -path [your_path] -output_name [your_output_name]```

The ```-path``` argument expects a folder with the following structure:
```
[your_path]
├── group 1
│   ├── image1.tif
│   ├── image2.tif
├── group 2
│   ├── image3.tif
│   ├── image4.tif
├── group 3
    ├── image5.tif
    ├── image6.tif
```

## Run segmentation_4.py
```python segmentation_4.py -path [your_path] -output_name [your_output_name]```

The ```-path``` argument expects a folder with the following structure:
```
[your_path]
├── folder1
│   ├── group 1
│   │   ├── image1.tif
│   │   ├── image2.tif
│   ├── group 2
│       ├── image3.tif
│       ├── image4.tif
└── folder2
    ├── group 3
    │   ├── image5.tif
    │   ├── image6.tif
    ├── group4
        ├── image7.tif
        ├── image8.tif
```
