// <copyright file="DataType.cs" company="Marek Behalek, marek.behalek@vsb.cz ,VSB-TU Ostrava">
// Copyright (c) Marek Behalek, marek.behalek@vsb.cz ,VSB-TU Ostrava. All rights reserved.
// </copyright>

namespace om38to13
{
    /// <summary>
    ///  Enum defining which data wil be used for testng.<br />
    /// </summary>
    public enum DataType
    {
        /// <summary>All data will be used.</summary>
        AllData,
        /// <summary>Only data from alignments will be used.</summary>
        JustAlignments,
        /// <summary> Only data from built assemblies will be used.</summary>
        JustAssemblies,
    }
}
