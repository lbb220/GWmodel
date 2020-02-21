$client = new-object System.Net.WebClient
$client.DownloadFile('https://haokunt-download-data.oss-cn-hangzhou.aliyuncs.com/GWmodelCUDA.zip', 'GWmodelCUDA.zip')
