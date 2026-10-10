(function exposeOpenQpHessian(globalScope) {
  "use strict";

  function firstNonEmptyArray(...values) {
    return values.find((value) => Array.isArray(value) && value.length)
      || values.find(Array.isArray)
      || null;
  }

  function numericValue(value) {
    if (value === null || value === undefined) return NaN;
    if (typeof value === "number") return Number.isFinite(value) ? value : NaN;
    if (typeof value !== "string") return NaN;
    const normalized = value.trim().replace(/[dD]/g, "E");
    if (!normalized) return NaN;
    if (/i$/i.test(normalized)) {
      const imaginary = Number(normalized.slice(0, -1));
      return Number.isFinite(imaginary) ? -Math.abs(imaginary) : NaN;
    }
    return Number(normalized);
  }

  function atomCountFromSource(source) {
    if (!source || typeof source !== "object") return 0;
    if (Array.isArray(source.atoms)) return source.atoms.length;
    if (Array.isArray(source.coord)) return Math.floor(source.coord.length / 3);
    for (const nested of [source.molecule, source.geometry, source.structure]) {
      const count = atomCountFromSource(nested);
      if (count) return count;
    }
    return 0;
  }

  function rawVectorComponents(rawMode) {
    if (!Array.isArray(rawMode) || !rawMode.length) return [];
    if (rawMode.every((vector) => vector && typeof vector === "object" && !Array.isArray(vector))) {
      return rawMode.flatMap((vector) => [vector.x, vector.y, vector.z].map(numericValue));
    }
    return rawMode.flat(Infinity).map(numericValue);
  }

  function normalizeModeVectors(rawMode, atomCount) {
    if (!atomCount) return [];
    const flat = rawVectorComponents(rawMode);
    if (flat.length < atomCount * 3) return [];
    const vectors = [];
    for (let index = 0; index < atomCount; index += 1) {
      const x = flat[index * 3];
      const y = flat[index * 3 + 1];
      const z = flat[index * 3 + 2];
      vectors.push({
        x: Number.isFinite(x) ? x : 0,
        y: Number.isFinite(y) ? y : 0,
        z: Number.isFinite(z) ? z : 0
      });
    }
    const maxLength = Math.max(...vectors.map((vector) => Math.hypot(vector.x, vector.y, vector.z)), 0);
    if (!maxLength) return [];
    return vectors.map((vector) => ({
      x: vector.x / maxLength,
      y: vector.y / maxLength,
      z: vector.z / maxLength
    }));
  }

  function emptyVibrationData() {
    return { modes: [], units: { frequency: "cm-1", ir: "km/mol", raman: "a.u." }, metadata: {} };
  }

  function modeRowsFromSource(source) {
    const rows = source?.vibrations?.modes;
    return Array.isArray(rows) && rows.some((row) => row && typeof row === "object" && !Array.isArray(row))
      ? rows
      : [];
  }

  function extractVibrations(source, explicitAtomCount) {
    if (!source || typeof source !== "object") {
      return emptyVibrationData();
    }
    const rows = modeRowsFromSource(source);
    const frequencies = firstNonEmptyArray(
      source.freqs,
      source.frequencies,
      source.frequency_modes?.["frequencies_cm-1"],
      source.vibrations?.frequencies,
      rows.length ? rows.map((row) => row.frequency ?? row.freq) : null
    );
    const rawModes = firstNonEmptyArray(
      source.modes,
      source.normal_modes,
      source.frequency_modes?.normal_mode_eigenvectors,
      rows.length ? rows.map((row) => row.vectors ?? row.vector ?? row.displacements ?? []) : source.vibrations?.modes
    );
    const infrared = firstNonEmptyArray(
      source.infrared_intensities,
      source.ir_intensities,
      source.vibrations?.ir,
      rows.length ? rows.map((row) => row.ir ?? row.infrared_intensity) : null
    );
    const raman = firstNonEmptyArray(
      source.raman_activities,
      source.vibrations?.raman,
      rows.length ? rows.map((row) => row.raman ?? row.raman_activity) : null
    );
    const atomCount = explicitAtomCount || atomCountFromSource(source);
    const modes = (frequencies || []).map((rawFrequency, index) => {
      const frequency = numericValue(rawFrequency);
      const ir = numericValue(infrared?.[index]);
      const ramanActivity = numericValue(raman?.[index]);
      return {
        index: Number(rows[index]?.index ?? rows[index]?.mode ?? index + 1),
        frequency,
        imaginary: frequency < 0,
        ir: Number.isFinite(ir) ? ir : null,
        raman: Number.isFinite(ramanActivity) ? ramanActivity : null,
        vectors: normalizeModeVectors(rawModes?.[index] ?? [], atomCount)
      };
    }).filter((mode) => Number.isFinite(mode.frequency));
    const metadata = source.vibrational_intensity_metadata || source.vibrations?.metadata || {};
    return {
      modes,
      units: {
        frequency: "cm-1",
        ir: metadata.ir_units || source.vibrations?.units?.ir || "km/mol",
        raman: metadata.raman_units || source.vibrations?.units?.raman || "a.u."
      },
      metadata
    };
  }

  function extractVibrationsFromLog(text, explicitAtomCount) {
    if (typeof text !== "string" || !text.trim()) return emptyVibrationData();
    const lines = text.split(/\r?\n/);
    let tableRows = [];

    lines.forEach((line, lineIndex) => {
      if (!/Mode\s+Frequency\s*\(\s*cm(?:\^?-?)1\s*\)/i.test(line)) return;
      const nextRows = [];
      for (let index = lineIndex + 1; index < lines.length; index += 1) {
        const row = lines[index].match(/^\s*(\d+)\s+([-+0-9.DEd]+i?)\s+([-+0-9.DEd]+)\s+([-+0-9.DEd]+)\s*$/i);
        if (!row) {
          if (nextRows.length) break;
          continue;
        }
        nextRows.push({
          index: Number(row[1]),
          frequency: numericValue(row[2]),
          ir: numericValue(row[3]),
          raman: numericValue(row[4])
        });
      }
      if (nextRows.length) tableRows = nextRows;
    });

    if (!tableRows.length) {
      lines.forEach((line) => {
        const row = line.match(/PyOQP\s+freq\s+(\d+)\s*:\s*(\S+)/i);
        if (!row) return;
        tableRows.push({
          index: Number(row[1]),
          frequency: numericValue(row[2]),
          ir: NaN,
          raman: NaN
        });
      });
    }

    let modeSectionStart = -1;
    lines.forEach((line, index) => {
      if (/Normal mode eigenvectors\s*\(Cartesian,\s*mass-unweighted\)/i.test(line)) {
        modeSectionStart = index;
      }
    });

    const vectorModes = [];
    if (modeSectionStart >= 0) {
      for (let index = modeSectionStart + 1; index < lines.length; index += 1) {
        const frequencyLine = lines[index].match(/^\s*Frequencies\s*--\s*(.+)$/i);
        if (!frequencyLine) continue;
        const frequencies = frequencyLine[1]
          .trim()
          .split(/\s+/)
          .map(numericValue)
          .filter(Number.isFinite);
        if (!frequencies.length) continue;

        const precedingTokens = (lines[index - 1] || "").trim().split(/\s+/);
        const modeIndexes = precedingTokens.length === frequencies.length
          && precedingTokens.every((token) => /^\d+$/.test(token))
          ? precedingTokens.map(Number)
          : frequencies.map((_frequency, offset) => vectorModes.length + offset + 1);

        let atomHeader = index + 1;
        while (atomHeader < Math.min(lines.length, index + 5) && !/\bAtom\s+AN\b/i.test(lines[atomHeader])) {
          atomHeader += 1;
        }
        if (atomHeader >= lines.length || !/\bAtom\s+AN\b/i.test(lines[atomHeader])) continue;

        const atomRows = [];
        let rowIndex = atomHeader + 1;
        for (; rowIndex < lines.length; rowIndex += 1) {
          const atomRow = lines[rowIndex].match(/^\s*(\d+)\s+(\d+(?:\.0+)?)\s+(?:[A-Za-z]{1,3}\s+)?(.+)$/);
          if (!atomRow) break;
          const components = atomRow[3]
            .trim()
            .split(/\s+/)
            .map(numericValue);
          if (components.length < frequencies.length * 3 || components.some((value) => !Number.isFinite(value))) break;
          atomRows.push(components);
        }

        if (atomRows.length) {
          frequencies.forEach((frequency, column) => {
            vectorModes.push({
              index: modeIndexes[column],
              frequency,
              rawVectors: atomRows.flatMap((components) => components.slice(column * 3, column * 3 + 3)),
              atomCount: atomRows.length
            });
          });
          index = rowIndex - 1;
        }
      }
    }

    const hasFrequencyTable = tableRows.length > 0;
    const rawModes = hasFrequencyTable ? tableRows : vectorModes;
    const modes = rawModes.map((row, position) => {
      const vectorMode = hasFrequencyTable
        ? vectorModes.find((mode) => mode.index === row.index)
        : vectorModes.find((mode) => mode.index === row.index) || vectorModes[position];
      const frequency = row.frequency;
      const ir = numericValue(row.ir);
      const raman = numericValue(row.raman);
      return {
        index: Number(row.index ?? vectorMode?.index ?? position + 1),
        frequency,
        imaginary: frequency < 0,
        ir: Number.isFinite(ir) ? ir : null,
        raman: Number.isFinite(raman) ? raman : null,
        vectors: normalizeModeVectors(
          vectorMode?.rawVectors || [],
          explicitAtomCount || vectorMode?.atomCount || 0
        )
      };
    }).filter((mode) => Number.isFinite(mode.frequency));

    return {
      modes,
      units: { frequency: "cm-1", ir: "km/mol", raman: "a.u." },
      metadata: modes.length ? { source: "OpenQP log" } : {}
    };
  }

  function extractHessianSummary(source) {
    if (!Array.isArray(source?.hessian)) return null;
    const dimension = source.hessian.length;
    let maxAbs = 0;
    let trace = 0;
    source.hessian.forEach((row, rowIndex) => {
      if (!Array.isArray(row)) return;
      row.forEach((value) => {
        const numeric = numericValue(value);
        if (Number.isFinite(numeric)) maxAbs = Math.max(maxAbs, Math.abs(numeric));
      });
      const diagonal = numericValue(row[rowIndex]);
      if (Number.isFinite(diagonal)) trace += diagonal;
    });
    return {
      dimension,
      maxAbs,
      trace,
      metadata: source.hessian_metadata || {}
    };
  }

  const api = {
    atomCountFromSource,
    extractHessianSummary,
    extractVibrations,
    extractVibrationsFromLog,
    normalizeModeVectors,
    numericValue
  };

  globalScope.OpenQPHessian = api;
  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  }
}(typeof globalThis !== "undefined" ? globalThis : this));
