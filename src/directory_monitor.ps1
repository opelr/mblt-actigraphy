## Monitor the Actigraphy CSV directory and run a nightly cron-job to update
## files changed or added in the last 24 hours

# Schtasks /create /tn "Actigraphy Data ETL" /sc daily /st 22:00 /tr "D:\lab\code\r\LabProjects\human\h4085\mblt-actigraphy\Scripts\actigraphy_csv_monitor.ps1"

cd D:\lab\code\r\LabProjects\human\h4085\mblt-actigraphy\
$Script_Path = "D:\lab\code\r\LabProjects\human\h4085\mblt-actigraphy\Scripts\00_actigraphy-cli.R"

$Data_Path = "D:\data\acquired_data\human\h4085\Actigraphy\CSV\"
$Filter = "*_Bedtime.csv"
$Yesterday = Get-Date (Get-Date).AddDays(-1) -UFormat "%Y-%m-%d"

$Updated_Files = Get-ChildItem -Path $Data_Path -Filter $Filter |
  Where-Object {$_.LastWriteTime -ge $Yesterday} |
  Get-ChildItem -Name

$Files_Join = $Updated_Files -join "---"
If ($Files_Join -eq "") {
  $Files_Join = " "
}

Rscript $Script_Path -f $Files_Join
