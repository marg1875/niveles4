# progreso.ps1 — Barras de progreso de la corrida de niveles4 (solo lectura).
# Uso:  powershell -ExecutionPolicy Bypass -File progreso.ps1         (una vez)
#       powershell -ExecutionPolicy Bypass -File progreso.ps1 -Vivo   (se refresca cada 30 s)
#       powershell -ExecutionPolicy Bypass -File progreso.ps1 -Vivo -Carpeta resultados_car
param([switch]$Vivo, [string]$Carpeta = 'resultados')
$raiz = Split-Path -Parent $MyInvocation.MyCommand.Path
$res  = Join-Path $raiz $Carpeta
$log  = Join-Path $res 'pipeline.log'

function Barra($nombre, $hecho, $total) {
    if ($total -le 0) { $total = 1 }
    $p = [math]::Min([double]1.0, [double]$hecho / [double]$total)
    $n = [int]($p * 30)
    "{0,-28} [{1}{2}] {3,3}%  ({4}/{5})" -f $nombre, ('#' * $n), ('.' * (30 - $n)), [int]($p*100), $hecho, $total
}

do {
    Clear-Host
    Write-Host "niveles4 - progreso  $(Get-Date -Format 'HH:mm:ss')`n"
    $ses = 0
    if (Test-Path "$res\sesiones.csv") { $ses = (Import-Csv "$res\sesiones.csv").Count }
    $sen = @(Get-ChildItem "$res\senales\*.mat" -ErrorAction SilentlyContinue | Where-Object Name -notlike '*.tmp.mat').Count
    Barra 'Etapa 1  preprocesamiento' $sen ([math]::Max($ses, 26))
    $car = @(Get-ChildItem "$res\caracteristicas\caracteristicas_d*.csv" -ErrorAction SilentlyContinue).Count
    Barra 'Etapa 2  caracteristicas' $car 9
    $arf = @(Get-ChildItem "$res\arff\d*\indice.csv" -ErrorAction SilentlyContinue).Count
    Barra 'Etapa 3  ARFF' $arf 9
    # Fase A: un archivo de cache por (duracion x dataset); 66 datasets por duracion
    $cacheA = @(Get-ChildItem "$res\weka\cache\d*\*.mat" -ErrorAction SilentlyContinue | Where-Object Name -notlike '*.tmp.mat').Count
    $totA = 0
    Get-ChildItem "$res\arff\d*\indice.csv" -ErrorAction SilentlyContinue | ForEach-Object { $totA += (Import-Csv $_.FullName).Count }
    if ($totA -eq 0) { $totA = 594 }
    Barra 'Etapa 4  Weka fase A' ([math]::Min($cacheA, $totA)) $totA
    $txt = if (Test-Path $log) { Get-Content $log -Raw -Encoding UTF8 } else { '' }
    $faseB = if ($txt -match 'fase B') { if ($txt -match 'Chequeo de fuga|Etapa 4: \d+ resultados') { 1 } else { 0.5 } } else { 0 }
    Barra 'Etapa 4  Weka fase B' ([int]($faseB*2)) 2
    $e5 = if ($txt -match 'Etapa 5: reportes') { 1 } else { 0 }
    Barra 'Etapa 5  reportes' $e5 1
    # Etapas 6/7: una barra por cada corrida de ajuste (resultados\ajuste*) con
    # su búsqueda (archivos *.avance.txt, 1120 modelos por sesión) y sus pasos.
    $fin = 1
    $dirs = @(Get-ChildItem "$res" -Directory -Filter 'ajuste*' -ErrorAction SilentlyContinue | Where-Object Name -ne 'ajuste_prueba')
    foreach ($d in $dirs) {
        $marca = 'conjunto .* -> ' + [regex]::Escape($d.FullName) + '\s*$'
        $ms = [regex]::Matches($txt, 'Etapa 6: ajuste[^\n]*' )
        $t6 = ''
        for ($k = $ms.Count - 1; $k -ge 0; $k--) { if ($ms[$k].Value -match $marca -or ($d.Name -eq 'ajuste' -and $ms[$k].Value -notmatch '->')) { $t6 = $txt.Substring($ms[$k].Index); break } }
        if ($t6.Length -gt 10) { $sig = $t6.IndexOf('Etapa 6: ajuste', 10); if ($sig -gt 0) { $t6 = $t6.Substring(0, $sig) } }
        $pasos = 0
        foreach ($m in @('diagn.stico listo', 'b.squeda lista', 'CV repetida lista', 'CV anidada de \w+ con', 'Etapa 6 lista')) { if ($t6 -match $m) { $pasos++ } }
        $hechos = 0; $total = 0
        $av = @(Get-ChildItem "$($d.FullName)\cache\busqueda_*.avance.txt" -ErrorAction SilentlyContinue)
        foreach ($a in $av) { $v = (Get-Content $a.FullName -TotalCount 1) -split ' '; $hechos += [int]$v[0]; $total += [int]$v[1] }
        if ($av.Count -lt 13) { $total += (13 - $av.Count) * 1120 }
        if ($av.Count -eq 0 -and $pasos -ge 2) { $hechos = 1; $total = 1 }   # corrida previa al contador
        if (Test-Path "$($d.FullName)\tabla_ajuste_kmi3c.csv") { $pasos = 5 }  # terminada (aunque el log se intercale)
        Barra ("{0} pasos" -f $d.Name) $pasos 5
        Barra ("{0} busqueda" -f $d.Name) $hechos $total
        if ($pasos -lt 5) { $fin = 0 }
    }
    # Etapa 8: reposo y eventos (pasos del log + datasets evaluados en la caché)
    $i8 = $txt.LastIndexOf('Etapa 8: reposo/eventos en 25')
    if ($i8 -lt 0) { $i8 = $txt.LastIndexOf('Etapa 8: reposo/eventos en 2') }
    if ($i8 -ge 0) {
        $t8 = $txt.Substring($i8); $p8 = 0
        foreach ($m in @('eventos y caracter.sticas listos', 'datasets ARFF', 'evaluaci.n lista', 'Etapa 8 lista')) { if ($t8 -match $m) { $p8++ } }
        Barra 'Etapa 8  reposo (pasos)' $p8 4
        $nd = 0; if ($t8 -match '(\d+) datasets ARFF') { $nd = [int]$Matches[1] }
        $nc = @(Get-ChildItem "$res\eventos\cache\*.mat" -ErrorAction SilentlyContinue | Where-Object Name -notlike '*.tmp.mat').Count
        if ($nd -gt 0) { Barra 'Etapa 8  datasets evaluados' $nc $nd }
        if ($p8 -lt 4) { $fin = 0 }
    }
    # Comparación de variantes (puentes x rechazo), en resultados_var\
    $var = Join-Path $raiz 'resultados_var'
    if (Test-Path $var) {
        foreach ($m in @('interpolar', 'eliminar')) {
            $ns = @(Get-ChildItem "$var\puentes_$m\senales\*.mat" -ErrorAction SilentlyContinue | Where-Object Name -notlike '*.tmp.mat').Count
            Barra ("Variantes: ICA $m") $ns 26
        }
        $nv = @(Get-ChildItem "$var\ventanas_*.mat" -ErrorAction SilentlyContinue).Count
        Barra 'Variantes: ventanas' $nv 3
        $c = if (Test-Path "$var\comparacion_variantes.csv") { 1 } else { 0 }
        Barra 'Variantes: 20 combinaciones' $c 1
        $lc = Join-Path $raiz 'resultados_car\pipeline.log'
        $tp = if (Test-Path $lc) { Get-Content $lc -Raw -Encoding UTF8 } else { '' }
        $np = ([regex]::Matches($tp, 'permutaci.n P\d-R\d')).Count
        Barra 'Variantes: permutaciones' ([math]::Min($np, 2)) 2
        if (-not (Test-Path "$var\permutacion_loso.csv")) { $fin = 0 }
        # Comparación ASR (comparar_asr.m)
        if (Test-Path "$var\asr_k20") {
            foreach ($k in @(20, 30)) {
                $ns = @(Get-ChildItem "$var\asr_k$k\senales\*.mat" -ErrorAction SilentlyContinue | Where-Object Name -notlike '*.tmp.mat').Count
                Barra ("ASR k=$k  preproceso") $ns 26
            }
            $nv = @(Get-ChildItem "$var\ventanas_asr*.mat" -ErrorAction SilentlyContinue).Count
            Barra 'ASR  ventanas' $nv 2
            $c = if (Test-Path "$var\comparacion_asr.csv") { 1 } else { 0 }
            Barra 'ASR  evaluacion variantes' $c 1
            $hp = 0
            foreach ($a in @(Get-ChildItem "$var\permutacion_*.avance.txt" -ErrorAction SilentlyContinue)) { $hp += [int]((Get-Content $a.FullName -TotalCount 1) -split ' ')[0] }
            Barra 'ASR  permutaciones (4x1000)' $hp 4000
            if (-not (Test-Path "$var\permutacion_final.csv")) { $fin = 0 }
        }
    }
    # Etapa 9: algoritmos fractales de Niveles3
    $f9 = Join-Path $res 'fractales'
    if (Test-Path $f9) {
        $leer = { param($n) $a = Join-Path $f9 "$n.avance.txt"; if (Test-Path $a) { (Get-Content $a -TotalCount 1) -split ' ' } else { @(0, 1) } }
        $p = & $leer 'pasos'
        Barra 'Etapa 9  pasos' ([int]$p[0]) ([int]$p[1])
        $na = @(Get-ChildItem "$res\arff\fractales\intra\*.arff" -ErrorAction SilentlyContinue).Count
        Barra 'Etapa 9  ARFF 3c intra' $na 13
        $vw = if (Test-Path "$f9\ventanas_fractales.mat") { 1 } else { 0 }
        Barra 'Etapa 9  ventanas contiguas' $vw 1
        $ng = @(Get-ChildItem "$f9\cache\gen_*.mat" -ErrorAction SilentlyContinue | Where-Object Name -notlike '*.tmp.mat').Count
        Barra 'Etapa 9  generalizacion' $ng 40
        $pp = & $leer 'permutacion'
        if ([int]$pp[0] -gt 0) { Barra 'Etapa 9  permutacion' ([int]$pp[0]) ([int]$pp[1]) }
        if ([int]$p[0] -lt [int]$p[1]) { $fin = 0 }
    }
    # Etapa 10: calibración de Higuchi y DFA
    $k10 = Join-Path $res 'calib_hig_dfa'
    if (Test-Path $k10) {
        $leer10 = { param($n) $a = Join-Path $k10 "$n.avance.txt"; if (Test-Path $a) { (Get-Content $a -TotalCount 1) -split ' ' } else { @(0, 1) } }
        $p = & $leer10 'pasos'
        Barra 'Etapa 10 pasos' ([int]$p[0]) ([int]$p[1])
        $na = @(Get-ChildItem "$res\arff\calib_hig_dfa\intra\*.arff" -ErrorAction SilentlyContinue).Count
        Barra 'Etapa 10 ARFF 3c intra' $na 13
        $nc = 0
        foreach ($c in @(Get-ChildItem "$k10\cache\c3_*.mat" -ErrorAction SilentlyContinue)) { $nc++ }
        Barra 'Etapa 10 3c (sesiones con cache)' $nc 13
        $l = & $leer10 'loso'
        Barra 'Etapa 10 LOSO (lote actual)' ([int]$l[0]) ([int]$l[1])
        $ng = @(Get-ChildItem "$k10\cache\gen_*.mat" -ErrorAction SilentlyContinue | Where-Object Name -notlike '*.tmp.mat').Count
        if ($ng -gt 0) { Barra 'Etapa 10 generalizacion' $ng 11 }
        foreach ($a in @(Get-ChildItem "$k10\perm_*.avance.txt" -ErrorAction SilentlyContinue)) {
            $v = (Get-Content $a.FullName -TotalCount 1) -split ' '
            Barra ('Etapa 10 ' + $a.BaseName.Replace('.avance', '').Substring(0, [math]::Min(22, $a.BaseName.Length - 7))) ([int]$v[0]) ([int]$v[1])
        }
        if ([int]$p[0] -lt [int]$p[1]) { $fin = 0 }
    }
    # Etapa 11: conjuntos ganadores de Niveles3
    $g11 = Join-Path $res 'niveles3_ganadores'
    if (Test-Path $g11) {
        $a = Join-Path $g11 'pasos.avance.txt'
        $p = if (Test-Path $a) { (Get-Content $a -TotalCount 1) -split ' ' } else { @(0, 5) }
        Barra 'Etapa 11 pasos' ([int]$p[0]) ([int]$p[1])
        $na = @(Get-ChildItem "$res\arff\niveles3_ganadores\intra3c\*.arff" -ErrorAction SilentlyContinue).Count
        Barra 'Etapa 11 ARFF combinados' $na 13
        $n3 = @(Get-ChildItem "$g11\cache\c3c_*.mat" -ErrorAction SilentlyContinue).Count
        $n2 = @(Get-ChildItem "$g11\cache\c2c_*.mat" -ErrorAction SilentlyContinue).Count
        Barra 'Etapa 11 intra 3c + 2c' ($n3 + $n2) 26
        $ng = @(Get-ChildItem "$g11\cache\gen_*.mat" -ErrorAction SilentlyContinue | Where-Object Name -notlike '*.tmp.mat').Count
        Barra 'Etapa 11 generalizacion' $ng 21
        foreach ($f in @(Get-ChildItem "$g11\perm_*.avance.txt" -ErrorAction SilentlyContinue)) {
            $v = (Get-Content $f.FullName -TotalCount 1) -split ' '
            Barra ('Etapa 11 permutacion ' + $f.BaseName.Split('_')[1].Split('.')[0] + ' (max)') ([int]$v[0]) ([int]$v[1])
        }
        if ([int]$p[0] -lt [int]$p[1]) { $fin = 0 }
    }
    # Etapa 12: laboratorio Higuchi/DFA
    $l12 = Join-Path $res 'laboratorio'
    if (Test-Path $l12) {
        $a = Join-Path $l12 'pasos.avance.txt'
        $p = if (Test-Path $a) { (Get-Content $a -TotalCount 1) -split ' ' } else { @(0, 9) }
        Barra 'Etapa 12 pasos' ([int]$p[0]) ([int]$p[1])
        foreach ($L in 250, 375, 512, 700, 1000) {
            $c3 = if (Test-Path "$l12\cache\A_W$L.mat") { 'listo' } elseif (Test-Path "$l12\cache\c3_W$L.mat") { 'curvas' } else { '-' }
            $f = Join-Path $l12 "bc_W$L.avance.txt"
            if (Test-Path $f) { $v = (Get-Content $f -TotalCount 1) -split ' '; Barra ("  W$L LOSO+exc (3c: $c3)") ([int]$v[0]) ([int]$v[1]) }
            elseif ($c3 -ne '-') { Barra ("  W$L LOSO+exc (3c: $c3)") 0 1 }
        }
        foreach ($f in @(Get-ChildItem "$l12\perm_*.avance.txt" -ErrorAction SilentlyContinue)) {
            $v = (Get-Content $f.FullName -TotalCount 1) -split ' '
            Barra ('  ' + $f.BaseName.Replace('.avance', '') + ' (max)') ([int]$v[0]) ([int]$v[1])
        }
        if ([int]$p[0] -lt [int]$p[1]) { $fin = 0 }
    }
    # Etapa 13: complejidad (todas las familias)
    $x13 = Join-Path $res 'complejidad'
    if (Test-Path $x13) {
        $a = Join-Path $x13 'pasos.avance.txt'
        $p = if (Test-Path $a) { (Get-Content $a -TotalCount 1) -split ' ' } else { @(0, 9) }
        Barra 'Etapa 13 pasos' ([int]$p[0]) ([int]$p[1])
        foreach ($L in 512, 700, 1000) {
            $c3 = Test-Path "$x13\cache\c3_W$L.mat"; $w = Test-Path "$x13\cache\w_W$L.mat"
            $nA = @(Get-ChildItem "$x13\cache\A_W${L}_*.mat" -ErrorAction SilentlyContinue | Where-Object Name -notlike '*.tmp.mat').Count
            if ($c3 -or $w -or $nA) {
                Barra ("  W$L curvas 3c / contiguas") ([int]$c3 + [int]$w) 2
                Barra ("  W$L fase 1 (3c, sesiones)") $nA 13
                $f = Join-Path $x13 "loso_W$L.avance.txt"
                if (Test-Path $f) { $v = (Get-Content $f -TotalCount 1) -split ' '; Barra ("  W$L fase 2 (LOSO)") ([int]$v[0]) ([int]$v[1]) }
                if (Test-Path "$x13\adelante_W$L.avance.txt") { Barra ("  W$L hacia adelante") 1 1 }
            }
        }
        $nf = @(Get-ChildItem "$x13\cache\final_*.mat" -ErrorAction SilentlyContinue).Count
        if ($nf -gt 0) { Barra '  validacion final (hechos)' $nf ([math]::Max($nf, 20)) }
        foreach ($f in @(Get-ChildItem "$x13\perm_*.avance.txt" -ErrorAction SilentlyContinue | Sort-Object LastWriteTime | Select-Object -Last 1)) {
            $v = (Get-Content $f.FullName -TotalCount 1) -split ' '
            Barra ('  ' + $f.BaseName.Replace('.avance', '') + ' (max)') ([int]$v[0]) ([int]$v[1])
        }
        if ([int]$p[0] -lt [int]$p[1]) { $fin = 0 }
    }
    # Etapa 15: controles de revisión
    $r15 = Join-Path $res 'revision'
    if (Test-Path $r15) {
        $a = Join-Path $r15 'pasos.avance.txt'
        $p = if (Test-Path $a) { (Get-Content $a -TotalCount 1) -split ' ' } else { @(0, 5) }
        Barra 'Etapa 15 pasos' ([int]$p[0]) ([int]$p[1])
        $nv = @(Get-ChildItem "$r15\cache\val_*.mat" -ErrorAction SilentlyContinue).Count
        Barra '  LOSO/LOPO/voto (configs)' $nv 9
        $no = @(Get-ChildItem "$r15\cache\orden_*.mat" -ErrorAction SilentlyContinue).Count
        Barra '  orden vs nivel (configs)' $no 9
        foreach ($f in @(Get-ChildItem "$r15\perm_*.avance.txt" -ErrorAction SilentlyContinue)) {
            $v = (Get-Content $f.FullName -TotalCount 1) -split ' '
            Barra ('  ' + $f.BaseName.Replace('.avance', '') + ' (max)') ([int]$v[0]) ([int]$v[1])
        }
        if ([int]$p[0] -lt [int]$p[1]) { $fin = 0 }
    }
    # Etapa 7: calibración
    $i7 = $txt.LastIndexOf('Etapa 7: calibraci')
    if ($i7 -ge 0) {
        $t7 = $txt.Substring($i7); $p7 = 0
        foreach ($m in @('ARFF de calibraci', 'candidatos individuales evaluados', 'pares evaluados', 'MEJOR:', 'CV anidada de la calibraci')) { if ($t7 -match $m) { $p7++ } }
        Barra 'Etapa 7  calibracion' $p7 5
        if ($p7 -lt 5) { $fin = 0 }
    }
    if (Test-Path $log) { Write-Host "`nUltimas lineas del log:"; Get-Content $log -Tail 4 -Encoding UTF8 }
    if ($Vivo -and -not $fin) { Start-Sleep -Seconds 30 }
} while ($Vivo -and -not $fin)
