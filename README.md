# watershed_segmenter
Scripts for instance segmentation of adipocytes 

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

## Run new_batch_segmentation.py
```python new_batch_segmentation.py -path [your_path] -output_name [your_output_name]```

The ```-path``` argument expects a ZIP folder with the following structure:
```
[your_path]
├── group 1
│   ├── image1.tif
│   ├── image2.tif
├── group 2
    ├── image3.tif
    ├── image4.tif
```

## Run segmentation_3.py
```python segmentation_3.py -path [your_path] -output_name [your_output_name]```

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
