; Inno Setup script for the NeoEdit Windows installer.
;   ISCC.exe /DVersion=0.2.0 packaging\neoedit.iss     (after pyinstaller packaging\neoedit.spec)
; Installs per-user by default (no administrator password needed), with an "install for all
; users" option; Start Menu entry, optional desktop icon and file associations, uninstaller.
#ifndef Version
  #define Version "0.0.0"
#endif
#define AppName "NeoEdit"
#define AppExe "NeoEdit.exe"

[Setup]
AppId={{7C1A2E5B-9F47-4D3A-8B1E-2E6F0A9C41D7}
AppName={#AppName}
AppVersion={#Version}
AppVerName={#AppName} {#Version}
AppPublisher=NeoEdit project
AppPublisherURL=https://github.com/evozoa/NeoEdit
AppSupportURL=https://github.com/evozoa/NeoEdit/issues
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir=..\dist
OutputBaseFilename=NeoEdit-Windows-Setup
SetupIconFile=..\src\neoedit\resources\icons\neoedit.ico
UninstallDisplayIcon={app}\{#AppExe}
UninstallDisplayName={#AppName}
LicenseFile=..\LICENSE
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ChangesAssociations=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"
Name: "assoc"; Description: "Open sequence files (.fasta, .bio, .gb, .aln, ...) with NeoEdit"; GroupDescription: "File associations:"

[Files]
Source: "..\dist\NeoEdit\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{group}\Uninstall {#AppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Registry]
Root: HKA; Subkey: "Software\Classes\NeoEdit.Sequence"; ValueType: string; ValueName: ""; ValueData: "Sequence alignment file"; Flags: uninsdeletekey; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\NeoEdit.Sequence\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\{#AppExe},0"; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\NeoEdit.Sequence\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\{#AppExe}"" ""%1"""; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.fasta"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.fasta\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.fas"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.fas\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.fa"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.fa\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.fna"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.fna\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.faa"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.faa\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.aln"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.aln\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.bio"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.bio\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.gb"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.gb\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.gbk"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.gbk\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.genbank"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.genbank\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.embl"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.embl\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.phy"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.phy\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.phylip"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.phylip\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.nex"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.nex\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.nexus"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.nexus\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.sto"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.sto\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.msf"; ValueType: string; ValueName: ""; ValueData: "NeoEdit.Sequence"; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.msf\OpenWithProgids"; ValueType: string; ValueName: "NeoEdit.Sequence"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc

[Run]
Filename: "{app}\{#AppExe}"; Description: "{cm:LaunchProgram,{#AppName}}"; Flags: nowait postinstall skipifsilent
