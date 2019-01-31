# CVWBM main readme

### Lots of materials about the project
1. Data sources table
2. Time Periods Table
3. Links to all the data
4. etc...


# Download NHD Data : 
```
mkdir nhd
cd nhd

curl -o 1802.zip https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_1802_HU4_Shape.zip && mkdir 1802 && tar -xvf 1802.zip && mv Shape/ 1802/ && rm 1802.zip

curl -o 1803.zip https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_1803_HU4_Shape.zip && mkdir 1803 && tar -xvf 1803.zip && mv Shape/ 1803/ && rm 1803.zip

curl -o 1804.zip https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_1804_HU4_Shape.zip && mkdir 1804 && tar -xvf 1804.zip && mv Shape/ 1804/ && rm 1804.zip

```